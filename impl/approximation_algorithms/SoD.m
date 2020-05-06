classdef SoD < AbstractApproximationAlgorithm   
    properties
        sub_size % function to determine that subset size
    end
    
    methods(Static)
        function val = getFileName()
            val = which(mfilename);
        end
                    
        function b = isDeterministic()
            b = false;
        end

    end
    
    
    methods  
        % Returns a name to be used for a file name.
        function name = getName(obj)
            name = 'SoD';
        end
        
        % Returns a name that can use latex formatting to be used in tikz
        % legend.
        function name = getLatexName(obj)
            name = 'SoD';
        end

    end
    
    methods (Access=protected)
        function initialize_impl(obj)  
            N = size(obj.X, 1);
            f = obj.K.get_number_of_inducing_inputs_f();
            obj.sub_size = @(step) ceil((N * (f(step))^2)^(1/3));
        end
        
        function [msg, t, times, t_pred, t_vfe, ymu_hat, ys2_hat, nlZ_hat, penalty, residual, sol_dist] = next_step_impl(obj, record_step)
            % we already gave the hyper-parameters to the kernel
            msg = 0;
            if ~record_step
                return
            end
            t = tic();
            cov     = @(varargin) obj.k(varargin{2:end});
            hyp.cov = zeros(size(obj.X, 2)+1, 1);
            hyp.lik = log(obj.sn2) / 2;
            m       = obj.sub_size(obj.current_step);
            try
                [ymu_hat, ys2_hat, nlZ_hat, ~] = gp(hyp, @infExact, [], cov, [], obj.X(1:m, :), obj.y(1:m), obj.Xs);
            catch ME
                warning(ME.message);
                msg = obj.SKS_CHOLESKY_FAILED;
            end
            t = toc(t);
            times = [];
            t_pred = 0;
            t_vfe = 0;
            penalty = 0;
            residual = NaN;
            sol_dist = NaN;
        end
    end
end

