function K = covRBF(hyp, x, varargin)

if nargin<2, K = '(2)'; return; end              % report number of parameters
K = covSEard([hyp(1) * ones(size(x, 2), 1); hyp(2)], x, varargin{:});