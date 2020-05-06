classdef AbstractApproximationAlgorithm < handle
    % Abstract class that defines which methods approximations algorithms
    % have to provide.
    properties (Constant = true)
        SKS_CHOLESKY_FAILED = 1;
        SKKS_CHOLESKY_FAILED = 2;
    end
    
    properties
        X % inputs
        y % targets
        K % kernel matrix
        trK % trace of K
        time_trK % wc-time to compute the trace of K
        k % kernel
        sn2 % noise
        Xs % test data points
        current_step; % the number of CG steps that have been executed so far
        steps % the number of CG steps that will be executed in total
        extended_recording; % whether to an extended recording or not
    end
    
    methods (Abstract, Static)
        getFileName() % returns the filename of the matlab file implement this class
        isDeterministic() % whether the algorithm is deterministic or not
    end  
    
    methods (Abstract)
        % Returns a name to be used for a file name.
        getName(obj)
        
        % Returns a name that can use latex formatting to be used in tikz
        % legend.
        getLatexName(obj)                
    end
    
    methods (Abstract, Access = protected)
        % initialization for the concrete algorithm
        initialize_impl(obj);
        
        % Executes the next step, record_step indiciates whether the
        % results of that step are recorded or not
        next_step_impl(obj, record_step)
    end
    
    methods
        function initialize(obj, X, y, K, k, sn2, Xs, steps, extended_recording)
            % As long as we do not change the objects below, no copy will
            % be performed.
            obj.X   = X;
            obj.y   = y;
            obj.K   = K;
            t = tic();
            obj.trK = sum(diag(K));
            obj.time_trK = toc(t);
            obj.k   = k;
            obj.sn2 = sn2;
            obj.Xs  = Xs;
            % initialize also serves as a reset. We must set the following values here.
            obj.current_step = 0;
            obj.steps = steps;
            obj.extended_recording = extended_recording;
            
            obj.initialize_impl();
        end
        
        function [msg, time, times, time_predict, time_vfe, ymu_hat, ys2_hat, nlZ_hat, penalty, residual, sol_dist] = next_step(obj, record_step)
            %
            % msg - 
            % t - accumulated time
            % times - additionaly algorithm dependent timings
            % 
            obj.current_step = obj.current_step + 1;
            if ~record_step
                msg = obj.next_step_impl(record_step);
                return;
            end
            
            if ~obj.extended_recording
                [msg, time, times, time_predict, time_vfe, ymu_hat, ys2_hat, nlZ_hat, penalty] = obj.next_step_impl(record_step);
                residual = -1;
                sol_dist = -1;
            else
                [msg, time, times, time_predict, time_vfe, ymu_hat, ys2_hat, nlZ_hat, penalty, residual, sol_dist] = obj.next_step_impl(record_step);
            end
        end
        
        function val = getHash(obj)
            % Takes some letters the document as hash value.
            
            content = fileread(obj.getFileName());
            % Remove white space.
            content(ismember(content, sprintf(' \t\n'))) = [];

            m       = numel(content);
            indices = [floor(m/4), ceil(m/4), floor(m/2), ceil(m/2), floor(3*m/4), ceil(3*m/4)];
            val     = content(indices);
        end
        
        function plot_overhead(obj, store_struct)
            % The purpose of this function is to visualize all the timing
            % measurements an algorithm takes.
            figure(1); clf; hold on;
            title('cumulative time for each step');
            xlabel('step');
            ylabel('time');

            %plot(store_struct.recorded_values(:, 5), store_struct.recorded_values(:, 1), 'DisplayName', '');
            plot(store_struct.steps, store_struct.recorded_values(4, :), 'DisplayName', 'self-measured training time');
            plot(store_struct.steps, store_struct.times(1, :), 'DisplayName', 'outside-loop time');
            plot(store_struct.steps, store_struct.recorded_values(4, :) + store_struct.recorded_values(5, :), 'DisplayName', 'train + prediction time');
            
            legend();
            
        end
    end
end

