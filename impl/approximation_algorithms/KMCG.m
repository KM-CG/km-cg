classdef KMCG < AbstractApproximationAlgorithm
    properties
        kXs; % k(X, Xs)
        dkXsXs; % diag(k(Xs, Xs))
        time_kernel; % The time need to evaluate above expressions.
        
        first_catch = Inf; % first time an exception is encountered
        
        cg_impl; % The CG implementation that will be used.
    end
        
    methods (Static)
        function val = getFileName()
            val = which(mfilename);
        end
        
        function b = isDeterministic()
            b = true;
        end
    end
    
    methods
        function obj = KMCG(cg_impl)
            obj.cg_impl = cg_impl;
        end
        
        function name = getName(obj)
            name = ['hybrid_'  obj.cg_impl.getName()];
        end
        
        function name = getLatexName(obj)
            name = ['\kmcg{}\cgImpl{' obj.cg_impl.getLatexName() '}'];
        end
    end
    
    methods(Access=protected)        
        function initialize_impl(obj)            
            t                = tic();
            obj.kXs          = obj.k(obj.X, obj.Xs);
            obj.dkXsXs       = obj.k(obj.Xs, 'diag');
            obj.time_kernel  = toc(t);
            
            obj.cg_impl.initialize(obj.X, obj.y, obj.K, obj.k, 0, obj.Xs, obj.steps, obj.extended_recording)
        end
        
        function [msg, t, times, t_pred, t_vfe, ymu_hat, ys2_hat, nlZ_hat, penalty, residual, sol_dist] = next_step_impl(obj, record_step)
            j = obj.current_step;
            N = numel(obj.y);
            
            if ~record_step
                msg = obj.execute_next_cg_step(j);
                return;
            end
            
            [msg, cg_cum_time, overhead_time, copy_time, S, Y, SY, cholSY] = obj.execute_next_cg_step(j);
            if msg ~= 0
                % CG crashed which makes it hard to continue.
                t        = [];
                times    = [];
                t_pred   = [];
                t_vfe    = [];
                ymu_hat  = [];
                ys2_hat  = [];
                nlZ_hat  = [];
                penalty  = [];
                residual = [];
                sol_dist = [];
                return
            end
                        
            t1 = tic(); % measure time
            Phib    = Y' * obj.y;
            logDetG = 2 * sum(log(diag(cholSY)));
            L       = Y' * Y + obj.sn2 * SY; L = 0.5 * L + 0.5 * L';
            try
                L           = chol(L);
                logDetSig   = 2 * sum(log(diag(L)));
                logDetSub   = logDetSig - logDetG + (N - j) * log(obj.sn2);
                fit         = L' \ Phib;
                fit         = (obj.y' * obj.y - fit' * fit) / obj.sn2;                
            catch ME
                warning(ME.message);

                msg = bitor(msg, AbstractApproximationAlgorithm.SKKS_CHOLESKY_FAILED);
                if j < obj.first_catch
                    obj.first_catch = j;
                    warning('SKKS stopped being spd after %d steps', j);
                end
                temp2       = cholSY' \ Y';
                logDetSub   = sum(log(diag(chol(temp2 * temp2' + obj.sn2 * eye(j))))) + (N - j) * log(obj.sn2);
                fit         = (obj.y' * obj.y - Phib' / L * Phib) / obj.sn2;
            end
            
                        
            nlZ_hat = fit + logDetSub + N * log(2 * pi);
            nlZ_hat = nlZ_hat / 2;
            
            t1 = toc(t1);
            
            % The overhead_time is not zero for text-book CG when computing
            % the Cholesky of S'Y. It's overhead for CG but NOT for this
            % algorithm.
            t     = cg_cum_time + overhead_time + t1;

            [t_pred, ymu_hat, ys2_hat, alpha_small] = obj.predict(j, cholSY, L, Phib, S);
            
            clear S; % to save memory
            
            % Titsias' penalty
            t_vfe   = tic();
            %temp2   = cholSY' \ Y'; this exceeds the memory limits
            penalty = 0.5 * (obj.trK - sum(sum((cholSY' \ Y').^2))) / obj.sn2;
            t_vfe   = toc(t_vfe) + obj.time_trK;
            
            times = [cg_cum_time; overhead_time; t1; copy_time];
            
            % No need to calculate the following quantities if we don't run
            % an extended recording.
            if ~obj.extended_recording, return; end
            
            alpha_large = (obj.y - Y * alpha_small) / obj.sn2;
            
            residual = norm(obj.y - obj.K * alpha_large - obj.sn2 * alpha_large);
            sol_dist = alpha_large' * (obj.K * alpha_large) + obj.sn2 * (alpha_large' * alpha_large) - 2 * (alpha_large' * obj.y);
        end
        
        function [msg, cg_cum_time, overhead_time, copy_time, S, Y, SY, cholSY] = execute_next_cg_step(obj, j)
            [msg, cg_cum_time, overhead_time, copy_time, S, Y, SY, cholSY] = obj.cg_impl.next_cg_step(j); 
        end
    end
    
    methods 
        % This function is not Access restricted so it can be tested
        function [t_pred, ymu_hat, ys2_hat, alpha_small] = predict(obj, j, cholSY, L, Phib, S)
            % predict
            t_pred = tic();
            SkXs   = S' * obj.kXs;
            if j < obj.first_catch
                alpha_small = solve_chol(L, Phib); 
                %ys2SoR      = L' \ SkXs;
                %ys2SoR      = sum(ys2SoR.^2, 1)' * obj.sn2;
                ys2SoR      = sum((L' \ SkXs).^2, 1)' * obj.sn2;
            else
                alpha_small = L \ Phib;
                ys2SoR      = obj.sn2 * diag(SkXs' / L * SkXs);
            end
            ymu_hat = SkXs' * alpha_small;
            ys2_hat = obj.dkXsXs - sum((cholSY' \ SkXs).^2, 1)' + ys2SoR + obj.sn2;
            
            t_pred  = toc(t_pred) + obj.time_kernel;
        end
    end
end

