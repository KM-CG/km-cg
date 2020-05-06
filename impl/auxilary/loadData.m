function [trainX, trainY, testX, testY] = loadData(name, fold, preprocess, folds, seed)
% loadData adapted from Krzysztof Chalupka http://homepages.inf.ed.ac.uk/ckiw/code/gpr_approx.html developed for:
% A Framework for Evaluating Approximation Methods for Gaussian Process Regression by Krzysztof Chalupka, Christopher K. I. Williams and Iain Murray (2012)
    EXPERIMENT.DATASET = name;
    EXPERIMENT.DATASET_FOLD = fold;
    EXPERIMENT.PREPROCESS_DATASET = preprocess;

    mydir = fullfile(pwd, 'data');

    restore_seed = randi(32000); % seed is set in the beginning and restored in the end
    if nargin < 5 
        seed = 0; %just that we always take the same shuffle
    end
    rng(seed);


    %rng('default');
    %----------------------------------------
    % MODIFY THIS TO
    % ADD YOUR OWN DATASET PREPROCESSING.
    %----------------------------------------
    if strcmp(EXPERIMENT.DATASET, 'MPG')
        mpg = csvread(fullfile(mydir, 'MPG/cleaned_auto-mpg.data'), 1);
        [n, D] = size(mpg);
        p = randperm(n); % seed is set in the beginning and restored in the end
        mpg = mpg(p, :);
        n_test = ceil(n / 10);
        n = n - n_test;
        trainX = mpg(1:n, 2:D);
        trainY = mpg(1:n, 1);
        testX = mpg(n+1:n+n_test, 2:D);
        testY = mpg(n+1:n+n_test, 1);
        clear mpg;
    elseif strcmp(EXPERIMENT.DATASET, 'CPU')
        cpuset = load(fullfile(mydir, 'CPU/cpu.mat'));
        trainX = cpuset.Xtrain';
        trainY = cpuset.ytrain';
        testX = cpuset.Xtest';
        testY = cpuset.ytest';
    elseif strcmp(EXPERIMENT.DATASET, 'PUMADYN')
        puma=load(fullfile(mydir, 'PUMADYN/pumadyn32nm.mat'));
        trainX = puma.X_tr;
        trainY = puma.T_tr;
        testX = puma.X_tst;
        testY = puma.T_tst;
        clear puma;
    elseif strcmp(EXPERIMENT.DATASET, 'PRECIPITATION')
        if folds == 0
            % this is the default setting when no cross-validation is desired
            EXPERIMENT.DATASET_FOLDS = 10;
            folds = 10;
        end
        prec = load(fullfile(mydir, 'PRECIPITATION/USprec1.txt'));
        stats = load(fullfile(mydir, 'PRECIPITATION/USprec2.txt'));
        trainY = sum(prec(prec(:,14)==0,2:13),2);
        trainY = trainY/100;
        trainX = stats(prec(:,14)==0,2:3);
        testX = []; % because folds is always greater 0 testX is initialized at the end
        testY = [];

        clear prec;
        clear stats;
    elseif strcmp(EXPERIMENT.DATASET, 'ABALONE')
        abalone = csvread(fullfile(mydir, 'ABALONE/abalone.data'), 1);
        N = size(abalone, 1);
        p = randperm(N);
        n = floor(N / 4 * 3);
        trainX = abalone(p(1:n), 1:end-1);
        trainY = abalone(p(1:n), end);
        testX = abalone(p(n+1:end), 1:end-1);
        testY = abalone(p(n+1:end), end);
        clear abalone;
    elseif strcmp(EXPERIMENT.DATASET, 'AILERONS')
        ailerons = csvread(fullfile(mydir, 'AILERONS/ailerons.data'), 1);
        trainX = ailerons(:, 1:end-1);
        trainY = ailerons(:, end);
        clear ailerons;
        ai_test = csvread(fullfile(mydir, 'AILERONS/ailerons.test'), 1);
        testX = ai_test(:, 1:end-1);
        testY = ai_test(:, end);
        clear ai_test;
    elseif strcmp(EXPERIMENT.DATASET, 'ELEVATORS')
        elevators = csvread(fullfile(mydir, 'ELEVATORS/elevators.data'), 1);
        trainX = elevators(:, 1:end-1);
        trainY = elevators(:, end);
        clear elevators;
        el_test = csvread(fullfile(mydir, 'ELEVATORS/elevators.test'), 1);
        testX = el_test(:, 1:end-1);
        testY = el_test(:, end);
        clear el_test;
     elseif strcmp(EXPERIMENT.DATASET, 'POLETELECOMM')
        pol = csvread(fullfile(mydir, 'POLETELECOMM/pol.data'), 1);
        trainX = pol(:, 1:end-1);
        trainY = pol(:, end);
        clear pol;
        pol_test = csvread(fullfile(mydir, 'POLETELECOMM/pol.test'), 1);
        testX = pol_test(:, 1:end-1);
        testY = pol_test(:, end);
        clear pol_test;   
    elseif strcmp(EXPERIMENT.DATASET, 'TOY')
        n = 400;
        % input locations
        X = [rand(n/4, 1); randn(n/4,1); 1 + 0.1*randn(n/4,1); -0.5 + 0.05 * randn(n/4,1)];
        trainX = sort(X(1:end/2));

        % kernel parameters
        sf2 = 1;        % amplitude
        ls  = .2;       % length scale
        sn2 = 0.02;     % noise

        % squared exponential kernel
        k   = @(a,b) sf2 * exp(- bsxfun(@minus,a,b').^2 ./ 2 ./ ls.^2);
        K = k(trainX,trainX);      % Gram matrix
        K = K + sn2 * eye(n/2);      % plus some noise
        L = chol(K);
        trainY = L' * randn(n/2,1);  % target for the linear program: sample from the GP.


        %res  = 300; % how many points to plot?
        %range= 3; % plot points from where to where?
        %testX  = 2 * range * ((1:res)' / res) - range;
        testX  = sort(X(end/2+1:end));
        testY  = (trainY' / L / L' * k(trainX, testX))';
        EXPERIMENT.PREPROCESS_DATASET = false;
    elseif strcmp(EXPERIMENT.DATASET, 'SOUND')
        if EXPERIMENT.PREPROCESS_DATASET
            EXPERIMENT.PREPROCESS_DATASET = false;
            warning('Preprocessing disabled for SOUND.');
        end
        if folds > 0
            folds = 0;
            warning('Cross validation disabled for SOUND as this would ruin the Toeplitz-property of the kernel matrix.');
        end
        
        testVars = {'xfull','xtrain','ytrain','xtest','ytest'};
        load('SOUND/audio_data.mat',testVars{:});
        trainX = xtrain;
        trainY = ytrain;
        testX = xtest;
        testY = ytest;
    elseif strcmp(EXPERIMENT.DATASET, 'TOEPLITZ_DEBUG_FULL') || strcmp(EXPERIMENT.DATASET, 'TOEPLITZ_DEBUG')
        trainX = [1:100, 151:200]';
        trainY = trainX;
        testX = (101:150)';
        testY = testX;
        if EXPERIMENT.PREPROCESS_DATASET
            EXPERIMENT.PREPROCESS_DATASET = false;
        end
    else
        error('Could not find dataset: %s!', name);
    end


    if folds > 1
        X = [trainX; testX];
        y = [trainY; testY];
        n = size(X, 1);
        p = randperm(n);
        X = X(p, :);
        y = y(p);

        a = (fold-1)*floor(n/folds) + 1;
        b = fold*floor(n/folds);
        testX = X(a:b, :);
        testY = y(a:b);
        trainX = X([1:(a-1),(b+1):n], :);
        trainY = y([1:(a-1),(b+1):n]);
    end

    rng('default');
    rng(restore_seed);

    if(~EXPERIMENT.PREPROCESS_DATASET)
        return;
    end
    %------------------------------------------------------------------------
    % DON'T MODIFY THIS.
    % Normalize data to zero mean, variance one. 
    %------------------------------------------------------------------------
    n = size(trainX, 1);
    trainYMean = mean(trainY);
    trainYStd  = std(trainY);
    %trainYStd(trainYStd == 0) = 1; % we don't want to divide by zero.
    stdTrainX = std(trainX);
    stdMatrix  = repmat(stdTrainX, n, 1);
    %stdMatrix(stdMatrix == 0) = 1;
    trainX = trainX(:, stdTrainX ~= 0); %remove useless features
    testX = testX(:, stdTrainX ~= 0);
    stdMatrix = stdMatrix(:, stdTrainX ~= 0);
    meanMatrix = repmat(mean(trainX), n, 1);

    D = size(trainX, 2); 
    trainX = (trainX - meanMatrix);
    trainX = trainX./stdMatrix;
    trainY = (trainY - trainYMean);
    trainY = trainY./trainYStd;

    testX  = (testX-repmat(meanMatrix(1,:), size(testX,1),1));
    testX = testX./repmat(stdMatrix(1,:), size(testX,1),1);
    testY  = (testY - trainYMean);
    testY = testY./trainYStd;

    if any(any(isnan(trainX) | isinf(trainX))), error('Training set contains NaN Values!'); end
    if any(any(isnan(trainY) | isinf(trainY))), error('Training targets contains NaN Values!'); end
    if any(any(isnan(testX) | isinf(testX))), error('Test set contains NaN Values!'); end
    if any(any(isnan(testY) | isinf(testY))), error('Training targets contains NaN Values!'); end
end
