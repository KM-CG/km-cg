% skeleton to make plots---can NOT be used on its own
% use show_one or save_*
DEBUG = false;

% define set of all algorithms ...
algos        = {FITC(), DTC(), KMCG(FOM()), FOM(), TextBookCG(), KMCG(TextBookCG()), KMCG(MatlabHookCG()), SoD(), KMCGmem(FOMmem()), DTC1Dgrid(), FITC1Dgrid()};
algos        = algos(algos_to_show); % ... and pick the ones we actually consider
colors       = [blu; blu; dre; mpg; mpg; dre; lightblu; cya; ora; gra'; 1 - lightblu; 1 - dre; 1 - blu; 1 - mpg; zeros(1, 3); 1 - cya; 1 - ora];
colors       = colors(algos_to_show, :);
lines        = {'--', ':','-', '-', '-.', '-.', '-.', '-.', '-.', '-.', '-.', '-.', '-.', '-.', '-.', '-.', '-.', '-.', '-.', '-.'};
lines        = lines(algos_to_show);
line_widths  = [1, 1, 1.4, 1.4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
line_widths  = line_widths(algos_to_show);

figure_hs    = cell(9,1); % cell to collect figure handles

if TIKZ
    close all;
    set(0, 'DefaultFigureVisible', 'off');
    data_sets_with_yaxis = {'ABALONE', 'MPG', 'ELEVATORS', 'linspace_distorted_100_100_', 'linspace_distorted_10_10_uniform_test_set'};
    if ~show_time
        data_sets_with_yaxis{end+1} = 'SOUND';
    end
    show_y_axis = ismember(data_set, data_sets_with_yaxis);
    
    data_sets_with_xaxis = {'SOUND', 'ELEVATORS', 'AILERONS', 'TOY'};
    show_x_axis = ismember(data_set, data_sets_with_xaxis);
else
    set(0, 'DefaultFigureVisible', 'on');
    show_x_axis = true;
    show_y_axis = true;
end


% get meta-data for the dataset
[~, ~, ~, X, ~, ~, ~, ~, ~, K] = loadOptimized(data_set, 0, false, 1, folds, cov_func);
[N, D]                         = size(X);
inducing_input_f               = K.get_number_of_inducing_inputs_f();
clear X K;

% common arguments for figure formating
format_figure_args = {TIKZ, show_time, show_x_axis, show_y_axis, data_set, N, D};
format_figure      = @(varargin) formater.format_figure(varargin{:});

figure_hs{6} = figure(6); clf; hold on; box on;% set(gcf, 'visible', ~TIKZ);
format_figure(format_figure_args{:});
if show_y_axis, ylabel('\smse{}'); end
title('smse');

figure_hs{7} = figure(7); clf; hold on; box on;% set(gcf, 'visible', ~TIKZ);
format_figure(format_figure_args{:});
if show_y_axis, ylabel('$\relerr{}$'); end
title('relative error');

figure_hs{8} = figure(8); clf; hold on; box on;% set(gcf, 'visible', ~TIKZ);
format_figure(format_figure_args{:});
if show_y_axis, ylabel('$\relerrVar{}$'); end
title('relative error of var');

figure_hs{9} = figure(9); clf; hold on; box on;% set(gcf, 'visible', ~TIKZ);
format_figure(format_figure_args{:});
if show_y_axis, ylabel('$\relerrNLZ{}$'); end
title('relative nlZ error');

% initialize cell of strings to store each algorithm's name
algo_names = cell([numel(algos) 1]);

max_time = 0;
min_time = Inf;

for fold = 1 : folds_to_show
for a = 1 : numel(algos)
    algo = algos{a};
    name      = algo.getName();
    algo_names{a} = algo.getLatexName();
    if strcmp(name, 'DTC')
        algo_names{a} = 'VFE'; % this is what we actually plot
    end
    
    
    % load all trials sequentially (as long as files exist)
    z = num2cell(zeros(1, 10));
    [perf, res, dist, msg, time, time_predict, rel, rel2, relZ, relZpenalty] = z{:};
    times          = [];
    recorded_steps = [];
    crashes        = [];
    rel_min   = Inf;
    rel_max   = 0;
    rel2_min  = Inf;
    rel2_max  = 0;
    relZ_min  = Inf;
    relZ_max  = 0;
    perf_min  = Inf;
    perf_max  = 0;
    vals_per_step = 0;
    n = 0;
    next_seed_exists = true;
    while next_seed_exists
        n    = n + 1;
        try
            loaded_file  = load(get_file_name(data_set, folds, fold, covfunc2str(cov_func), name, seed + n - 1));
            store_struct = loaded_file.store_struct;
            if ~isequal(algo.getHash(), store_struct.hash)
                if strcmp(data_set, 'SOUND') && (strcmp('FITC', algo.getName()) || strcmp('DTC', algo.getName()))
                    % in this case the results are from FITC1d or DTC1d
                else
                    warning('The loaded result file was created with another version of the algorithm (%s)! Action?', name);
                    keyboard
                end
            end
            
            s     = size(store_struct.recorded_values, 2);
            if s == 0, continue; end % can't do much if the algorithm crashes immediately
            crash = max_step + 1;
            if store_struct.msg ~= 0
                %s = s - 1;
                crash = s;
            end

            if s > length(recorded_steps), recorded_steps = store_struct.steps(1:s); end
            
            if recorded_steps(s) > max_step
                idx = find(recorded_steps <= max_step);
                s = length(idx);
                recorded_steps = store_struct.steps(1:s);
            end
            
            if crash <= s
                crashes = [crashes crash]; %#ok<AGROW>
            end
            
            recorded_values = store_struct.recorded_values(:, 1:s);
            
            t = recorded_values(4, :); % training times for the current seed
            
            
            nlZ_penalty = 0;
            if strcmp(name, 'DTC')
                % make DTC to VFE
                nlZ_penalty = recorded_values(9, :);
                t = t + recorded_values(end-1, :); % this is the time it took to compute the penalty term
            end
            
            vals_per_step = broadcast(@(a, b) a + b, vals_per_step, ones(1, s));
            v1 = vals_per_step(1:length(time));
            v2 = vals_per_step(1:s);
            % calculate the mean
            delta = broadcast(@(a, b) b - a, time ./ v1, t ./ v2);
            time  = broadcast(@(a, b) a + b, time, delta(1:length(v2)));
            times(n, 1:s) = t; %#ok<SAGROW>
            
            rel_n = recorded_values(6, :);           
            [rel, rel_min, rel_max] = compute_statistics(rel_n, rel, rel_min, rel_max, v1, v2);
           
            perf_n= recorded_values(1, :);
            [perf, perf_min, perf_max] = compute_statistics(perf_n, perf, perf_min, perf_max, v1, v2);
            
            rel2_n  = recorded_values(7, :);
            if strcmp(name, 'FOM')
                rel2_n = NaN(size(rel2_n));
            end
            [rel2, rel2_min, rel2_max] = compute_statistics(rel2_n, rel2, rel2_min, rel2_max, v1, v2);

            nlZ  = store_struct.nlZ;

            relZ_n = abs((recorded_values(8, :) - nlZ_penalty - nlZ) / nlZ);
            [relZ, relZ_min, relZ_max] = compute_statistics(relZ_n, relZ, relZ_min, relZ_max, v1, v2);           
        catch MatExc
            if n > 1 && strcmp(MatExc.identifier, 'MATLAB:load:couldNotReadFile')
                next_seed_exists = false;
            else
                rethrow(MatExc);
            end
        end
    end
    
    if ~show_time
        x_axis = recorded_steps;
    else
        times(times == 0) = NaN; %#ok<SAGROW>
        if n > 2
            med = median(times, 'omitnan');
        else
            med = times;
        end
        if strcmp(name, 'DTC') || strcmp(name, 'FITC')
            % fit quadratic function for algorithms that do not measure
            % their time sequentially
            phi    = @(x) [ones(size(x)) x x.^2];
            phix   = phi(inducing_input_f(recorded_steps)');
            x_axis = (phix * (phix'*phix \ (phix'*med')))';
        else
            x_axis = med;
        end
    end
    
    y_clip_value = 0;
    x_clip_value = Inf;
    if show_time
        if TIKZ
            y_clip_value = 1e-3;
        end
        
        if strcmp(name, 'FITC')
            max_time = x_axis(end);
        end
        min_time = min(min_time, min(med));
        
        x_clip_value = max_time;        
    end
    
    args_in = {lines{a}, colors(a, :), line_widths(a), x_clip_value, y_clip_value, crashes};
    plot_figure(figure_hs{6}, x_axis, perf, perf_min, perf_max, args_in{:});
    plot_figure(figure_hs{7}, x_axis, rel, rel_min, rel_max, args_in{:});
    plot_figure(figure_hs{8}, x_axis, rel2, rel2_min, rel2_max, args_in{:});
    plot_figure(figure_hs{9}, x_axis, relZ, relZ_min, relZ_max, args_in{:});
end
end

if show_time && max_time > 0
    % cut-off x-axis after FITC
    for j = 6 : 9
        set(0, 'currentfigure', figure_hs{j}); %figure(j);
        % exclude 0 because pgfplots always uses it as a tick
        xlim([min_time max_time]);
    end
end

figure_hs{4} = figure(4); clf; hold on; box off; axis off; 
algos_in_legend = numel(algos_to_show):-1:1;
format_legend;
title('legend');
