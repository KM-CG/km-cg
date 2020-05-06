%% create plots that demonstrate how unstable the textbook version of CG is
%% set parameters
seed         = 12345
TIKZ         = true;
data_sets = {'ABALONE', 'AILERONS', 'ELEVATORS', 'MPG', 'POLETELECOMM', 'PUMADYN', 'PRECIPITATION', 'TOY'};
cov_funcs = {@(varargin) covSEard(varargin{:})};

folds         = 2;
folds_to_show = 1;
algos_to_show = (3:6)';
max_step      = 100;
formater      = OnlyCGFormater();

prefix        = 'cg_unstable_';
show_time     = false;

%% iterate over the standard regression datasets
for d = 1 : numel(data_sets)
    data_set       = data_sets{d};
    for c = 1 : numel(cov_funcs)
        cov_func   = cov_funcs{c};
        show
        % we will save only the relative error w.r.t. the mean, that's
        % sufficient
        savetikz(figure_hs(7), [prefix data_set sprintf('_folds_%i_folds_visible_%i_', folds, folds_to_show) covfunc2str(cov_func)], seed );
        close all;
    end 
end

%% save legend
fh = figure(4); clf; hold on; box off; axis off;
algos_in_legend = numel(algos_to_show):-1:1;
format_legend;
title('legend');
savelegend(fh, prefix);
