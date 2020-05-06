%%
seed         = 12345
TIKZ         = true;
data_sets = {'linspace_distorted_10_10_uniform_test_set', 'linspace_distorted_1000_1000_uniform_test_set', 'linspace_distorted_100_100_100_uniform_test_set', 'linspace_distorted_30_30_30_30_uniform_test_set'};
cov_funcs = {@(varargin) covSEard(varargin{:})};
folds         = 2;
folds_to_show = 1;
algos_to_show = 1:4;
max_step      = 100;
formater      = StructuredMatrixFigureFormater();
format_figure = @(varargin) formater.format_figure(varargin{:});



prefix        = 'cg_hybrid_grid_';
show_time     = false;

%%
for d = 1 : numel(data_sets)
    data_set       = data_sets{d};
    
    for c = 1 : numel(cov_funcs)
        cov_func   = cov_funcs{c};
        show
        savetikz( figure_hs(6:9), [prefix data_set sprintf('_folds_%i_folds_visible_%i_', folds, folds_to_show) covfunc2str(cov_func)], seed );
        close all;
    end 
end

%%
prefix        = 'cg_hybrid_grid_time_';
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