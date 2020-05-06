% This script is a sceleton for the other run_experiment_* files.
for d = 1 : numel(data_sets)    
    for c = 1 : numel(cov_funcs)
        data_set       = data_sets{d};
        cov            = cov_funcs{c};
        [hyp, gp_time, smse, X, y, Xs, ys, ys2, nlZ, K, yKy] = loadOptimized(data_set{:}, fold, folds, cov);
        cf             = @(varargin) cov(hyp.cov, varargin{:}); 
        sn2            = exp(2 * hyp.lik);
        [N, D]         = size(X);   
        rng(seed); % set seed for each dataset and algorithm

        steps = min(5000, N); % do at most 5000 steps
        % if the dataset is larger than 100.000 datapoints restrict the number of steps to 500
        if N > 10^5, steps = 500; end 
        
        for a = 1 : numel(approx_algos)
            for s = 0 : (repetitions-1)
                rng(seed+s); % set seed for each dataset and algorithm
                % permuting the dataset is not too costly - it's okay to have it here in the loop
                p = randperm(N);
                X = X(p, :);
                y = y(p);
                %K = K(p, p);
                K.setPermutation(p);
                
                algo = approx_algos{a};
                % the result storage gets cov and NOT cf because it needs the
                % name
                resultStorage = ResultStorage(gp_time, algo.getName(), algo.getHash(), ys, ys2, nlZ, data_set{1}, folds, fold, cov, hyp, smse, yKy, seed+s);
                algo.initialize(X, y, K, cf, sn2, Xs, steps, extended_recording);
                time = 0;
                for step = 1 : steps
                    record_step = step < 200 || mod(step, 25) == 0; % record the first 200 steps and then every 25
                    if record_step
                        % give some debug output
                        fprintf('dataset %i, cov %i, algo %i, seed %i, step %i\n', d, c, a, s+1, step);
                        t = tic();
                        [msg, time, times, time_predict, time_vfe, ymu_hat, ys2_hat, nlZ_hat, penalty, residual, sol_dist] = algo.next_step(record_step);
                        t = toc(t);
                        resultStorage.store(step, ymu_hat, ys2_hat, nlZ_hat, penalty, residual, sol_dist, msg, time, [t; times(:)], time_predict, time_vfe);                        
                    else
                        msg = algo.next_step(record_step);
                    end
                    % Time is measured within the algorithm to allow some
                    % parts to be excluded from the measurements.
                    

                    if msg ~= 0
                        % We can not continue if the algorithm crashed.
                        break 
                    end 
                    
                    if time > 120 && step > 100, break; end % break when it becomes too expensive
                end
                    
                resultStorage.write();
                
                % reverse permutation of the dataset
                ip = zeros(N, 1);
                ip(p) = 1 : N;
                X = X(ip, :);
                y = y(ip);
                clear ip;
                %K = K(ip, ip); % no longer necessary because the
                %permutation is set in the beginning                
            end
        end
    end 
end
    