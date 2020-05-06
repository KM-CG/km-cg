classdef FITC < AbstractApproximationAlgorithm
    properties 
        num_ind % function that determines the number of inducing inputs
        dkXs
        kernel_time
        
        first_catch;
        second_catch;
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
            name = 'FITC';
        end
        
        % Returns a name that can use latex formatting to be used in tikz
        % legend.
        function name = getLatexName(obj)
            name = 'FITC';
        end
        
    end
    
    methods(Access=protected)
        function initialize_impl(obj)
            t = tic();
            obj.dkXs = obj.k(obj.Xs, 'diag');
            obj.kernel_time = toc(t);
            
            obj.num_ind = obj.K.get_number_of_inducing_inputs_f();
            
            % initialize also serves as a reset. We must set these values here.
            obj.first_catch = inf;
            obj.second_catch = inf;
        end
        
        function [KXU, KUU] = get_inducing_input_matrices(obj, m)
            KXU = obj.K(:, 1:m);
            KUU = KXU(1:m, :);
        end
        
        function KXSU = get_test_inducing_cov(obj, m)
            KXSU = obj.k(obj.X(1:m, :), obj.Xs);
        end
        
        function [msg, t, times, t_pred, t_vfe, ymu_hat, ys2_hat, nlZ_hat, penalty, residual, sol_dist] = next_step_impl(obj, record_step)
            msg = 0;
            if ~record_step, return; end
            t = tic();
            j = obj.current_step;
            N = numel(obj.y);
            m = obj.num_ind(j);
            [Phi, cholSY] = obj.get_inducing_input_matrices(m);
            % cholSY will be chol(K(1:m, 1:m)) later
            L = obj.sn2 * cholSY; % we want L to contain sn2 * K(1:m, 1:m)
            try
                cholSY      = chol(cholSY);                    
                logDetG     = 2 * sum(log(diag(cholSY)));

                dQff  = Phi / cholSY; dQff = sum(dQff.^2, 2);
            catch
                msg = AbstractApproximationAlgorithm.SKS_CHOLESKY_FAILED;
                if j < obj.first_catch
                    obj.first_catch = j;
                    warning('SKS stopped being spd after %d steps', j);
                end
                if j < 100 && N < 10000
                    % the dataset is small enough s.t. we can try to
                    % compute the solution the inefficient way
                    dQff     = diag(Phi / cholSY * Phi');
                    logDetG  = NaN;
                else
                    [t, times, t_pred, t_vfe, ymu_hat, ys2_hat, nlZ_hat, penalty, residual, sol_dist] = deal(nan);
                    return;
                end
            end
            sqrtLambda= sqrt((obj.k(obj.X, 'diag') - dQff) / obj.sn2 + 1);

            %Phi   = diag(1 ./ sqrt(Lambda)) * Phi;
            Phi   = bsxfun(@rdivide, Phi, sqrtLambda);
            b     = obj.y ./ sqrtLambda;
            L     = L + Phi' * Phi;
            Phib  = Phi' * b;


            try
                L           = chol(L);
                logDetSig   = 2 * sum(log(diag(L)));
                logDetSub   = logDetSig - logDetG + 2 * sum(log(sqrtLambda)) + (N - m) * log(obj.sn2);
                fit         = L' \ Phib;
                fit         = (b' * b - fit' * fit) / obj.sn2;

                alpha_small =  solve_chol(L, Phib); 

            catch
                if j < obj.second_catch
                    obj.second_catch = j;
                    warning('SKKS stopped being spd after %d steps', j);
                end
                msg = bitor(msg, AbstractApproximationAlgorithm.SKKS_CHOLESKY_FAILED);
                
                if j < 100 && N < 10000
                    % the dataset is small enough s.t. we can try to
                    % compute the solution the inefficient way
                    Km   = obj.K(:, 1:m);
	                Khat = Km / Km(1:m, :) * Km'; Khat = Khat / 2 + Khat' / 2;
	                Khat = Khat + (diag(obj.k(obj.X, 'diag') - diag(Khat))) + obj.sn2 * eye(N);
	
	                try
	                    Khat         = chol(Khat);
	                    logDetSub    = 2 * sum(log(diag(Khat)));
                    catch
	                    % well, we're screwed
	                    logDetSub    = NaN;
                    end
	                alpha_small  = L \ Phib;
	                fit          = Khat' \ b; fit = fit' * fit;
                else
                    [t, times, t_pred, t_vfe, ymu_hat, ys2_hat, nlZ_hat, penalty, residual, sol_dist] = deal(nan);
                    return;
                end
            end


            nlZ_hat   = fit + logDetSub + N * log(2 * pi);
            nlZ_hat   = nlZ_hat / 2;
            t         = toc(t);
            
            times     = [];
            
            t_pred = tic();        
            SkXs   = obj.get_test_inducing_cov(m);
            
            if j < obj.first_catch
                temp     = cholSY' \ SkXs;
                ys2prior = sum(temp.^2, 1)';
            else
                % since the first cholesky threw an exception, cholSY == K(1:m, 1:m)
                ys2prior = diag(SkXs' / cholSY * SkXs);
            end
            
            if j < obj.second_catch
                ys2SoR      = L' \ SkXs;
                ys2SoR      = sum(ys2SoR.^2, 1)' * obj.sn2;
            else
                ys2SoR       = obj.sn2 * diag(SkXs' / L * SkXs);   
            end
            

            ymu_hat   = SkXs' * alpha_small;
            ys2_hat   = obj.dkXs - ys2prior + ys2SoR + obj.sn2;
            
            t_pred    = toc(t_pred) + obj.kernel_time;
            
            % Titsias' penalty
            penalty = 0; % by construction
            t_vfe   = 0;
            
            if ~obj.extended_recording, return; end

            alpha_large = (obj.y - Phi * alpha_small) / obj.sn2;
            residual = norm(obj.y - obj.K * alpha_large - obj.sn2 * alpha_large);
            sol_dist = alpha_large' * (obj.K * alpha_large) + obj.sn2 * (alpha_large' * alpha_large) - 2 * (alpha_large' * obj.y);
        end
    end
end

