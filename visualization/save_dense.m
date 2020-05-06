startup;
%% set parameters
seed         = 12345
TIKZ         = true;
data_sets = {'ABALONE', 'AILERONS', 'ELEVATORS', 'MPG', 'POLETELECOMM', 'PUMADYN', 'PRECIPITATION', 'TOY'}; %
cov_funcs = {@(varargin) covMaternard(5, varargin{:}), @(varargin) covSEard(varargin{:})};

folds         = 2;
folds_to_show = 1;
algos_to_show = 1:4; % show VFE, FITC, CG and KMCG
max_step      = 100; % only the first 100 steps are actually interesting

formater      = DenseMatrixFigureFormater();
format_figure = @(varargin) formater.format_figure(varargin{:});


%% save figures where x-axis is the number of CG-steps
prefix        = 'cg_hybrid_comparison_';
show_time     = false;

for d = 1 : numel(data_sets)
    data_set       = data_sets{d};
    for c = 1 : numel(cov_funcs)
        cov_func   = cov_funcs{c};
        show % create figures
        % save smse, relative error of mean, of variance and approximation
        % error to log-marginal likelihood
        savetikz( figure_hs(6:9), [prefix data_set sprintf('_folds_%i_folds_visible_%i_', folds, folds_to_show) covfunc2str(cov_func)], seed );
        close all;
    end 
end

%% save figures plotted over training time
prefix        = 'cg_hybrid_comparison_time_';
show_time     = true;
max_step      = 5000;

for d = 1 : numel(data_sets)
    data_set       = data_sets{d};
    for c = 1 : numel(cov_funcs)
        cov_func   = cov_funcs{c};
        show
        savetikz( figure_hs(6:9), [prefix data_set sprintf('_folds_%i_folds_visible_%i_', folds, folds_to_show) covfunc2str(cov_func)], seed );
        close all;
    end 
end

%% save legends
fh = figure(4); clf; hold on; box off; axis off;
algos_in_legend = numel(algos_to_show):-1:1; %#ok<NASGU>
format_legend;
title('legend');
savelegend(fh, 'all'); % legend containing all algorithms

fh = figure(12); clf; hold on; box off; axis off;
algos_in_legend = [1, 2, 3];
format_legend;
title('legend_without_cg');
savelegend(fh, 'without_cg'); % legend without CG
