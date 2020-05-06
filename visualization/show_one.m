% show results for a single dataset and covariance function
seed         = 12345
TIKZ         = false;
data_set     = 'ABALONE';
folds        = 2;
folds_to_show= 1;
algos_to_show= 3:6; % which algorithms to show
cov_func     = @(varargin) covSEard(varargin{:}); %@(varargin) covMaternard(5, varargin{:});

formater      = DenseMatrixFigureFormater();
format_figure = @(varargin) formater.format_figure(varargin{:});

show_time     = false;
if ~show_time
    max_step = 100;
else
    max_step = 5000;
end


show

if TIKZ
    prefix = 'test_';
    savetikz( figure_hs([7, 8]), [prefix data_set sprintf('_folds_%i_folds_visible_%i_', folds, folds_to_show) covfunc2str(cov_func)], seed, true );
end