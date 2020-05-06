classdef ResultStorage < handle
    properties %(Access = protected)
        file_path_name
        msg = 0; % the last error message
        recorded_values
        times % different time-measurements per step (number of values may vary with each algorithm)
        num_recorded_values = 9; % the number of values that we record
        ys % test targets
        ys2 % corresponding variances
        nlZ % evidence
        absys % absolute values, where too small values have been replaced
        absys2 % absolute values of the posterior variance
        GPsmse % smse of the exact GP
        yKy % y / K * y
        t % the time the exact GP takes
        data_set % dataset name
        cov_func % name of the covariance function
        hyps % hyper-parameters
        seed % seed
        ysvar % variance of the test targets
        nys % norm of the test targets
        name % name of the algorithm
        hash % hash created from the algorithm's source file
        recorded_steps % The steps that get recorded
    end
    
    methods
        function obj = ResultStorage(t, name, hash, ys, ys2, nlZ, data_set, folds, fold, cov_func, hyps, smse, yKy, seed)
            obj.t               = t;
            obj.name            = name;
            obj.hash            = hash;
            obj.ysvar           = var(ys);
            obj.ys              = ys;
            obj.ys2             = ys2;
            obj.nlZ             = nlZ;
            obj.absys           = abs(ys); obj.absys(obj.absys < 1e-8) = 1e-8;
            obj.absys2          = abs(ys2); obj.absys2(obj.absys2 < 1e-8) = 1e-8;
            obj.nys             = norm(ys);
            obj.hyps            = hyps;
            obj.cov_func        = covfunc2str(cov_func);
            obj.GPsmse          = smse;
            obj.yKy             = yKy;
            obj.recorded_steps  = [];
            obj.recorded_values = [];
            obj.times           = [];
            obj.file_path_name  = get_file_name(data_set, folds, fold, obj.cov_func, name, seed);
        end
        
        function ex = check_if_results_exist(obj)
            try 
                load(obj.file_path_name, 'store_struct')
                if isfield(store_struct, 'hash') && isequal(obj.hash, store_struct.hash)
                    ex = true;
                    return
                else
                   % delete results and start again
                   warning('Deleting obsolete result file: %s.', obj.file_path_name)
                   delete([obj.file_path_name '.mat']);
                   %error('starting experiments from scratch');
                end
            catch ME
                warning(ME.message);
            end
            ex = false;
        end
        
        function obj = store(obj, step, mu, var, nlZ, penalty, residual, solDist, msg, time, times, time_predict, time_vfe)
            % This purpose of this function is to calculate and to store
            % all relevent quantities an approximaton algorithm may
            % produce.
            %
            % args
            % obj  - this object
            % step - The number of steps CG has taken.
            % mu1  - Mean prediction of the test targets.
            % var1 - Uncertainty over the test targets.
            % nlZ   - Negative marginal log-likelihood.
            % penalty - additional term derived by Titsias: Variational Learning of Inducing Variables in Sparse Gaussian Processes (2009)
            % residual - y - (K + sigma^2 I) * x
            % solDist - is the distance to the solution in A-norm.
            % msg - a number with a meaning defined by the algorithm 
            % time - time used for training as measured by the algorithm
            % times - additional measurements that are algorithm dependent, the first entry is the time measured outside of the loop
            % time_predict - time used for prediction
            % time_vfe - time used to compute the penalty term            
            if msg ~= 0
                obj.msg = msg;
                return
            end
            to_store                     = [obj.smse(mu), residual, solDist + obj.yKy, time, time_predict, obj.rel_err1(mu), ...
                obj.rel_err2(var), nlZ, penalty, 0, time_vfe];
            obj.recorded_values = [obj.recorded_values, to_store'];
            obj.times = [obj.times, times'];
            obj.recorded_steps = [obj.recorded_steps, step];
        end
        
        function write(obj)
            store_struct.msg             = obj.msg;
            store_struct.t               = obj.t;
            store_struct.nlZ             = obj.nlZ;
            store_struct.recorded_values = obj.recorded_values;
            store_struct.times           = obj.times;
            store_struct.steps           = obj.recorded_steps;
            store_struct.hyps            = obj.hyps;
            store_struct.hash            = obj.hash;
            store_struct.smse            = obj.GPsmse;
            store_struct.ysvar           = obj.ysvar; %#ok<STRNU>
            save(obj.file_path_name, 'store_struct');
        end
    end
    
    methods(Access = protected)
        function e = msll(obj, mu_hat, var_hat)
            e = 0.5 * mean( (obj.ys - mu_hat).^2 ./ var_hat ) + 0.5 * mean( log(2 * pi * var_hat) );
        end
        
        function e = smse(obj, mu_hat)
            e = mean( (obj.ys - mu_hat).^2 ) / obj.ysvar;
        end
        
        function e = rel_err1(obj, mu_hat)
            e = abs(obj.ys - mu_hat);
            e = mean(e ./ obj.absys);
        end
        
        function e = rel_err2(obj, mu_hat)
            % computes the average relative error for the posterior
            % variance
            e = abs(obj.ys2 - mu_hat);
            e = mean(e ./ obj.absys2);
        end
    end
end