classdef LanczosSR < AbstractCG
    % Lanczos with selected reorthogonalization
    properties
        alpha
        beta
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
        function name = getName(obj)
            name = 'LanczosSR';
        end
        
        function name = getLatexName(obj)
            name = '\lanczosSR{}';
        end
        
        function [msg, cum_time, overhead_time, copy_time, S, Y, SY, LSY] = next_cg_step(obj, j)
                t   = tic; % measure time
                msg = 0;
                SY  = obj.LanczosSRstep(j);
                if nargout > 2
                    overhead_time = tic;
                    S = obj.S(:, 1:j);
                    Y = obj.Y(:, 1:j);
                    copy_time = toc(overhead_time);
                    SY = S' * Y; SY = 0.5 * SY + 0.5 * SY';
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
                
                t = toc(t);
                
                obj.cum_time = obj.cum_time + t - overhead_time;
                cum_time     = obj.cum_time;
        end
        
        function initialize_cg(obj)
            t     = tic();
            N     = size(obj.X, 1);
            obj.S = zeros(N, obj.steps);
            obj.Y = zeros(N, obj.steps);

            obj.x  = zeros(N, 1);
            obj.r  = obj.y; % we start with x0=0
            obj.nr = norm(obj.r);
            obj.s  = obj.r / obj.nr; % normalize by norm
            obj.cum_time = toc(t);
        end
        
        function Ti = LanczosSRstep(obj, i)
            % Matlab code LanczosSelectOrthog.m
            % For "Applied Numerical Linear Algebra",  Chapter 7, section 5
            % Written by James Demmel, Jun  6, 1997
            %
            % Perform Lanczos with selective orthogonalization
            %
            % https://people.eecs.berkeley.edu/~demmel/ma221_Fall09/Matlab/LANCZOS_README.html
                  z = obj.K * obj.s + obj.sn2 * obj.s;
                  obj.S(:, i)  = obj.s;
                  obj.Y(:, i) = z;
                  obj.alpha(i) = obj.s'*z;
                  z = z - obj.alpha(i)*obj.s;
                  if i>1
                      z = z - obj.beta(i-1)*obj.S(:,i-1); 
                  end
                  obj.beta(i) = norm(z);
                  obj.s = z/obj.beta(i);
            %
            %     selectively orthogonalize, when the error bound < sqrt(eps)
                  if (i>1)
                     Ti = diag(obj.alpha(1:i)) + diag(obj.beta(1:i-1),1) + diag(obj.beta(1:i-1),-1);
                  else
                     Ti = obj.alpha(1);
                  end
                  [V,D] = eig(Ti); % Simon: supposedly Matlab is aware that the matrix is tridiagonal
            %     sort eigenvalues into decreasing order
                  [Ds,Is] = sort(-diag(D));
                  Vs = V(:,Is);
                  ErrorBoundsFLSO = abs((obj.beta(i)*Vs(i,:))');
                  select = find( ErrorBoundsFLSO <= sqrt(eps)*max(abs(Ds)) );
                  if ~isempty(select)
            %         for s = select'
            %           ritzvector = QFLSO*Vs(:,s);
            %           z = z - (ritzvector'*z)*ritzvector;
            %         end
                    ritzvector = obj.S(:,1:i)*Vs(:,select);
                    z = z - sum(bsxfun(@times, ritzvector, (ritzvector'*z)'), 2);
                    obj.beta(i) = norm(z);
                    obj.s = z/obj.beta(i);
                  end
            %
                  obj.S(:,i+1) = obj.s;
                  %RitzComponentsFLSO(1:i,i) = ((QFLSO(:,i+1)'*QFLSO(:,1:i))*Vs)';
        end
    end    
end

