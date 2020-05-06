classdef MatlabHookCG < AbstractCG
    properties
        cg_times % times that it costs to execute one step
        hook_j % internal step counter
        tic_stamp; % stamp to measure time between matrix multiplcations
        preconditioner = [];
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
        function obj = MatlabHookCG(varargin)
            obj = obj@AbstractCG(varargin{:});
        end
        
        function name = getName(obj)
            name = 'MatlabCG';
        end
        
        function name = getLatexName(obj)
            name = 'CG';
        end
        
        function setPreconditioner(obj, preconditioner)
            obj.preconditioner = preconditioner;
        end
        
        function [msg, cum_time, overhead_time, copy_time, S, Y, SY, LSY] = next_cg_step(obj, j)
                msg = 0;
                
                % this is not part of the vanilla CG implementation
                if nargout > 2
                    overhead_time = tic;
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
                cum_time = sum(obj.cg_times(1:j));
        end
        
        function initialize_cg(obj)
            N     = size(obj.X, 1);
            obj.S = zeros(N, obj.steps);
            obj.Y = zeros(N, obj.steps);
            obj.cg_times = zeros(1, obj.steps);
            % obj.cum_time = 0; we do not need this variable
            
            tol        = 0; % use matlab's default tolerance
            obj.hook_j = 0;
            
            obj.tic_stamp = tic;
            % execute CG using this object as hook (see function below)
            if isempty(obj.preconditioner)
                pcg(@(s) obj.mtimes(s), obj.y, tol, obj.steps);
            else
                pcg(@(s) obj.mtimes(s), obj.y, tol, obj.steps, obj.preconditioner);
            end
            obj.cg_times(obj.hook_j) = toc(obj.tic_stamp);
        end
        
        function y = mtimes(obj,s)
            if obj.hook_j == 0
                % in the very first step CG multiplies with zeros
                y = zeros(size(s));
                obj.hook_j = 1;
                return
            end
            y                    = obj.K * s + obj.sn2 * s;
            obj.S(:, obj.hook_j) = s;
            obj.Y(:, obj.hook_j) = y;
            obj.cg_times(obj.hook_j) = toc(obj.tic_stamp);
            obj.tic_stamp        = tic;
            obj.hook_j           = obj.hook_j + 1;
        end
    end
    
    
    methods(Access=protected)    
        function [msg, t, times, t_pred, t_vfe, t_test_nlZ, ymu_hat, ys2_hat, nlZ_hat, penalty, test_nlZ, residual, sol_dist] = next_step_impl(obj, record_step)
            [msg, t, overhead_time, copy_time, S, Y, SY, LSY] = obj.next_cg_step(obj.current_step);
            if ~record_step, return; end
            if msg
                [t, times, t_pred, t_vfe, t_test_nlZ, ymu_hat, ys2_hat, nlZ_hat, penalty, test_nlZ, residual, sol_dist] = deal(nan);
                return
            end
            obj.x    = S * solve_chol(LSY, S' * obj.y);
            nlZ_hat  = NaN;
            
            t_pred   = tic();
            ymu_hat  = obj.kXs' * obj.x;
            ys2_hat  = NaN(size(ymu_hat));
            t_pred   = toc(t_pred) + obj.time_kernel;
            
            times      = []; % no need for any extra information
            penalty    = NaN;
            t_vfe      = 0;
            t_test_nlZ = 0;
            test_nlZ   = NaN;
            
            if ~obj.extended_recording, return; end
            obj.nr   = norm(obj.y - obj.K * obj.x);
            residual = obj.nr; %norm(obj.r);
            sol_dist = obj.x' * (obj.K * obj.x) + obj.sn2 * (obj.x' * obj.x) - 2 * obj.x' * obj.y;
        end
    end    
end