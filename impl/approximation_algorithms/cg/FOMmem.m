classdef FOMmem < AbstractCG
    % memory efficient version of FOM
    methods (Static)
        function val = getFileName()
            val = which(mfilename);
        end
        
        function b = isDeterministic()
            b = true;
        end
    end
    
    methods
        function name = getName(obj)
            name = 'FOMmem';
        end
        
        function name = getLatexName(obj)
            name = '\FOM{}';
        end
    
        
        function [msg, cum_time, overhead_time, copy_time, S, Y, SY, LSY, alphas, betas] = next_cg_step(obj, j)
                t = tic; % measure time
                % sn2 = 0 for KMCG
                y           = obj.K * obj.s + obj.sn2 * obj.s;
                obj.S(:, j) = obj.s;%normalization is done in the try block
                
                copy_time = tic;
                S = obj.S(:, 1:j);
                copy_time = toc(copy_time);
                
                try
                    obj.SY(1:j, j) = S' * y;
                    obj.SY(j, 1:j) = obj.SY(1:j, j);
                    obj.SY(j, j)   = y' * obj.s;
                    SY    = obj.SY(1:j, 1:j);
                    LSY   = chol(SY);
                    alphas     = solve_chol(LSY, S' * obj.y);
                    obj.x = S * alphas;
                    obj.r = obj.y - obj.K * obj.x; % if mvm is cheap this is actually fast enough
                    obj.nr= norm(obj.r);
                    betas = -solve_chol(LSY, S' * (obj.K * obj.r));
                    obj.s = (obj.r + S * betas) / obj.nr;
                    msg   = 0;
                catch ME
                    warning(ME.message);
                    warning('Stopping FOM after %d steps instead of %d with residual norm %d.', j, obj.steps, obj.nr);
                    msg = AbstractApproximationAlgorithm.SKS_CHOLESKY_FAILED;
                    [cum_time, overhead_time, copy_time, S, Y, SY, LSY] = deal(nan);
                    return
                end
                t = toc(t);
                obj.cum_time = obj.cum_time + t;
                cum_time = obj.cum_time;
                overhead_time = 0;
                Y = [];
        end
        
        function initialize_cg(obj)
            t = tic();
            N     = size(obj.X, 1);
            obj.S = zeros(N, obj.steps);
            obj.SY= zeros(obj.steps, obj.steps);
            obj.r = obj.y;
            obj.nr= norm(obj.r);
            obj.s = obj.r / obj.nr;
            obj.cum_time = toc(t);
        end
    end 
end

