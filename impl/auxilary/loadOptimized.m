function [hyp, gp_time, smse, X, y, Xs, ys, ys2, nlZ, K, yKy] = loadOptimized(dataSetName, optimizationSteps, targetsAreExactGP, fold, folds, cov)
% loads a dataset and kernel hyper-parameters
    if startsWith(dataSetName, 'pv__step')
        % This is a prefix for some extra 'linspace' experiments.
        % Entering this if-branch means that the dataset is loaded for
        % visualization.
        underscores = strfind(dataSetName, '_');
        dataSetName = dataSetName(underscores(4)+1:end); 
    end

    if startsWith(dataSetName, 'linspace_distorted')
        % grid datasets get a special treatment
        [hyp, gp_time, smse, X, y, Xs, ys, ys2, nlZ, K, yKy] = load_grid_data_set(dataSetName(numel('linspace_distorted'):end), cov);
        return
    end

    preprocess = true;
    [X, y, Xs, ys] = loadData(dataSetName, fold, preprocess, folds);
    hyp.lik = log(1e-3) / 2;             % initial noise value before hyper-parameter tuning
    hyp.cov = zeros([size(X, 2)+1, 1]);  % length scales and amplitude
    gp_time = 0;                         % the time it took for computing the final Cholesky
    smse    = 0;                         % final standardized mean squared error of the GP for the dataset
    
    K = [];                              % initialize empty kernel matrix
    
    if strcmp(dataSetName, 'SOUND')
        if optimizationSteps > 0
            warning('Setting optimizationSteps=0 for SOUND.'); 
            optimizationSteps = 0;
        end
        if targetsAreExactGP
            warning('SOUND is too large for exact GP approximation - disabling.'); 
            targetsAreExactGP = false;
        end
        
        load('data/SOUND/audio_data.mat', 'hypc');
        hyp = hypc; % hyper-parameters are part of the dataset
        
        cf = @(varargin) cov(hyp.cov, varargin{:});
        k  = cf((X(1):X(end))', X(1));
        K  = Toeplitz(k, Xs); % initialize a Toeplitz matrix
    end
    
    if optimizationSteps > 0
        % tune hyper-parameters if desired
        hyp     = hyper_parameters(dataSetName, folds, fold, cov, X, y, hyp, optimizationSteps);
    end

    if targetsAreExactGP
        % predict with the exact GP and use the predictions as targets
        sn2     = exp(2 * hyp.lik);
        
        disp('Initializing full GP');
        gp_time = tic();
        [ymu, ys2, nlZ, post] = gp(hyp, @infExact, [], cov, [], X, y, Xs);
        gp_time = toc(gp_time);
        smse = sum((ymu - ys).^2) / size(ys, 1)
        if smse > 0.05, warning('The GP seems to have bad hyper-parameters!'); end
        ys = ymu; % Targets are the test predictions of the GP.
        ys2 = diag(ys2); % keep only the variances
        disp('Done.');
        
        clear post.K;
        
        yKy = post.L' \ (y / sn2); yKy = yKy' * yKy;
    else
        % if we do not compute the exact posterior there is no posterior
        % variance or marginal log-likelihood that we could use as
        % ground-truth
        ys2 = zeros(size(ys));
        yKy = 0;
        nlZ = 0;
    end
    
    if isempty(K)
        % if the kernel matrix has not been initialized yet, build a dense
        % matrix
        cf = @(varargin) cov(hyp.cov, varargin{:}); 
        K  = DenseMatrix(cf(X));
    end
end

function [hyp, gp_time, smse, X, y, Xt, ys, ys2, nlZ, K, yKy]  = load_grid_data_set(dataSetName, k )
% this function parses the parameters from the dataset name. The actual
% loading happens in the function below.
    gp_time = NaN;
    smse    = NaN;
    axis_points = @(d, g) linspace(-g/4, g/4, g)' + 1e-3 * randn([g, 1]);
    
    uniform_test_set = false;
    if endsWith(dataSetName, 'uniform_test_set')
        uniform_test_set = true;
        % remove the postfix s.t. it does not interfere with parsing the
        % dataset dimensions in the next part
        dataSetName = dataSetName(1:end-numel('uniform_test_set'));
    end
    
    % parse the dimensionality of each axis
    underscores = strfind(dataSetName, '_');
    tensor_size = [];
    for j = 1 : (numel(underscores) - 1)
        tensor_size = [tensor_size, str2double(dataSetName(underscores(j)+1:underscores(j+1)-1))]; %#ok<AGROW>
    end
    
    % use default hyper-parameters of length-scale 1 and noise 0.001
    cov_hyp     = zeros(size(tensor_size, 2) + 1, 1);
    sn2         = 1e-3;
    hyp.cov     = cov_hyp;
    hyp.lik     = log(sn2) / 2;
    
    if uniform_test_set
        Xt_fun         = @(X) (rand(100, numel(tensor_size)) - 0.5) * diag(tensor_size) / 2;    
        missing_points = [];
    else
        missing_points = randperm(prod(tensor_size)); 
        missing_points = missing_points(1:floor(end/2));
        Xt_fun         = @(X) X(missing_points, :);
    end
    [X, y, Xt, ys, ys2, nlZ, K, yKy] = load_grid_data_set_(tensor_size, axis_points, missing_points, Xt_fun, k, cov_hyp, sn2);
end

function [X, y, Xt, ys, ys2, nlZ, K, yKy] = load_grid_data_set_(tensor_size, axis_points, missing_points, Xt_fun, k, cov_hyp, sn2)
    D = numel(tensor_size);
    Np= prod(tensor_size);
    
    id = setdiff(1:Np, missing_points); % indices of points that will be part of the dataset

    % build input grid, Kronecker kernel matrix and compute Eigenvalues and
    % Eigenvectors to sample targets
    Ks    = cell([D, 1]);
    Xg    = cell(size(Ks)); % all inputs in form of a cell array
    EVecs = cell(size(Ks));
    evals = 1;
    for d = 1 : D
        g         = tensor_size(d);
        xd        = axis_points(d, g);
        Xg{d}     = xd;
        
        Ks{d}      = k([cov_hyp(d) cov_hyp(D+1)/D], xd);
        [Q, E]     = eig(Ks{d});
        EVecs{d}   = Q;
        evals      = evals * diag(E)'; % * evals';
        evals      = evals(:);
    end
    
    
    
    [y, ~] = get_sample(Np, EVecs, evals, sn2); % draw a sample from the full grid
    y      = y(id);
    
    
    [Xg{:}]= ndgrid(Xg{:});
    X      = zeros(Np, D); % will be resized
    for d = 1 : D
        X(:, d) = reshape(Xg{d}, [Np 1]);
    end
    
    % Xt_fun needs access to the full grid!
    Xt     = Xt_fun(X);
    
    % remove missing points from data-set
    X      = X(id, :); 
    
    
    K       = Kron(@(varargin) k(cov_hyp, varargin{:}), X, Ks, missing_points); % build Kronecker kernel matrix
    
    [ys, ys2, yKy, nlZ, ~, ~] = exact_GP(Np, X, y, EVecs, evals, k, cov_hyp, sn2, Xt); % compute exact posterior of the full grid
    assert(var(ys) > 1e-6); % let's make sure the ys are not all 0
    assert(all(ys2 > 0));
    assert(isreal(nlZ));
end

function [y, sample] = get_sample(Np, EVecs, evals, sn2)
    sample  = randn([Np 1]);
    
    % for the exact GP we also use the missing entries (otherwise it's not
    % possible).
    y       = MultiKronMVM(EVecs, Np, 1:Np, MultiKronMVMT(EVecs, Np, 1:Np, sample(:)) .* sqrt(evals + sn2));

end

function [yt, ys2, yKy, nlZ, kXt, solvSqrt] = exact_GP(Np, X, y, EVecs, evals, k, cov_hyp, sn2, Xt)
    N  = size(X, 1);
    Nt = size(Xt, 1);
    kXt= k(cov_hyp, X, Xt);
    if Np > N
        % some points are missing
        % we have to resort to the naive calculation...
        L      = chol(k(cov_hyp, X) + sn2 * eye(size(X, 1)));
        alpha  = L' \ y;
        temp   = L' \ kXt;
        logdet = 2 * sum(log(diag(L)));

        solvSqrt = @(varargin) error('solvSqrt is a debug output and not meant for actual use.');
    else
        solvSqrt= @(v) MultiKronMVM(EVecs, Np, 1:Np, MultiKronMVMT(EVecs, Np, 1:Np, v) ./ sqrt(evals + sn2));
        
        temp = zeros([N Nt]);
        for s = 1 : Nt
            temp(:, s) = solvSqrt(kXt(:, s));
        end
        alpha= solvSqrt(y);

        logdet = sum(log(evals + sn2));
    end
    yt   = temp' * alpha;
    ys2  = k(cov_hyp, Xt, 'diag') - sum(temp.^2)' + sn2;
    yKy  = alpha' * alpha;
    nlZ  = yKy / 2 + logdet / 2 + N * log(2 * pi) / 2;
end
