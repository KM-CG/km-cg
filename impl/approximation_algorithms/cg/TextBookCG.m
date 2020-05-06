classdef TextBookCG < AbstractCG
    % Basically the book version of CG and using     
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
            name = 'CG';
        end
        
        function name = getLatexName(obj)
            name = 'CG';
        end
        
        function [msg, cum_time, overhead_time, copy_time, S, Y, SY, LSY] = next_cg_step(obj, j)
                
                msg = 0;
                t   = tic; % measure time
                y   = obj.K * obj.s + obj.sn2 * obj.s;
                
                % this is not part of the vanilla CG implementation
                if nargout > 2
                    overhead_time = tic;
                    %nr  = norm(obj.r);
                    obj.Y(:, j) = y; % already normalized: % / nr; % normalize with norm of residual
                    obj.S(:, j) = obj.s; % / nr;
                    S = obj.S(:, 1:j);
                    Y = obj.Y(:, 1:j);
                    copy_time = toc(overhead_time);
                    SY  = S' * Y; SY = 0.5 * SY + 0.5 * SY';
                    try
                        LSY = chol(SY);
                    catch ME
                        warning(ME.message);
                        warning('Stopping CG after %d steps instead of %d with residual norm %d.', j, obj.steps, obj.nr);
                        msg = AbstractApproximationAlgorithm.SKS_CHOLESKY_FAILED;
                        LSY = [];
                    end
                    overhead_time = toc(overhead_time);
                else
                    overhead_time = 0;
                end
                                
                ys      = y' * obj.s;
                rs      = obj.r' * obj.s;
                alpha   = rs / ys;
                obj.x   = obj.x + alpha * obj.s;                
                obj.r   = obj.r - alpha * y;
                
                nr      = norm(obj.r);
                %beta    = obj.r' * y / ys;
                %obj.s   = obj.s / nr;
                obj.s   = obj.r / nr + nr / obj.nr * obj.s;
                obj.nr  = nr;
                t = toc(t);
                obj.cum_time = obj.cum_time + t - overhead_time;
                cum_time = obj.cum_time;
        end
        
        function initialize_cg(obj)
            N     = size(obj.X, 1);
            obj.S = zeros(N, obj.steps);
            obj.Y = zeros(N, obj.steps);

            t = tic();
            obj.tol = 0.01 * norm(obj.y);
            obj.x = zeros(N, 1);
            obj.r = obj.y; % we start with x0=0
            obj.nr = norm(obj.r);
            obj.s = obj.r / obj.nr; % normalize by norm
            obj.cum_time = toc(t);
        end
    end    
end

