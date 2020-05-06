seed         = 12345
TIKZ         = true;
data_set     = 'SOUND';

max_step      = 100;
formater      = StructuredMatrixFigureFormater();
folds         = 2;
folds_to_show = 1;
algos_to_show = 1:4;
    
cov_func   = @(varargin) covSEard(varargin{:});

show_time     = false;
prefix        = 'cg_hybrid_comparison_';
show
savetikz( figure_hs(6:7), [prefix data_set sprintf('_folds_%i_folds_visible_%i_', folds, folds_to_show) covfunc2str(cov_func)], seed );

close all;
max_step      = 5000;
show_time     = true;
prefix        = 'cg_hybrid_comparison_time_';
show
savetikz( figure_hs(6:7), [prefix data_set sprintf('_folds_%i_folds_visible_%i_', folds, folds_to_show) covfunc2str(cov_func)], seed );
