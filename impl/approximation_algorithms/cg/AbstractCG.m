classdef AbstractCG < AbstractApproximationAlgorithm    
    properties
        tol % allowed tolerance when to stop CG
        S
        Y
        SY % S'Y
        cum_time
        x
        s
        r
        nr % norm of r
        
        kXs;
        time_kernel; % The time need to evaluate above expression.

    end
    
    methods(Abstract)
        % Initializes the concrete CG algorithm. The rationale behind this
        % extra method is to distinguish the cases when used as GP
        % approximation algorithm or inside KMCG.
        initialize_cg(obj);

        
        % The step is part of the argument such that implementations can be
        % used without calls to #next_step().
        [msg, cum_time, overhead_time, copy_time, S, Y, SY, cholSY] = next_cg_step(obj, step);
    end
    
    methods(Access = protected)
        function initialize_impl(obj)    
            t                = tic();
            obj.kXs          = obj.k(obj.X, obj.Xs);
            obj.time_kernel  = toc(t);
            obj.initialize_cg();
        end
        
        function [msg, t, times, t_pred, t_vfe, ymu_hat, ys2_hat, nlZ_hat, penalty, residual, sol_dist] = next_step_impl(obj, record_step)
            [msg, t] = obj.next_cg_step(obj.current_step);
            if ~record_step, return; end
            nlZ_hat  = NaN;
            
            t_pred   = tic();
            ymu_hat  = obj.kXs' * obj.x;
            ys2_hat  = NaN(size(ymu_hat));
            t_pred   = toc(t_pred) + obj.time_kernel;
            
            times      = []; % no need for any extra information
            penalty    = NaN;
            t_vfe      = 0;
            
            if ~obj.extended_recording, return; end
            residual = obj.nr; %norm(obj.r);
            sol_dist = obj.x' * (obj.K * obj.x) + obj.sn2 * (obj.x' * obj.x) - 2 * obj.x' * obj.y;
        end
        
    end
    
    
end

