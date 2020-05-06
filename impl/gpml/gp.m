function [varargout] = gp(hyp, inf, mean, cov, lik, x, y, xs)
% Gaussian Process inference and prediction. The gp function provides a
% flexible framework for Bayesian inference and prediction with Gaussian
% processes for scalar targets, i.e. both regression and binary
% classification. The prior is Gaussian process, defined through specification
% of its mean and covariance function. The likelihood function is also
% specified. Both the prior and the likelihood may have hyperparameters
% associated with them.
%
% Two modes are possible: training or prediction: if no test cases are
% supplied, then the negative log marginal likelihood and its partial
% derivatives w.r.t. the hyperparameters is computed; this mode is used to fit
% the hyperparameters. If test cases are given, then the test set predictive
% probabilities are returned. Usage:
%
%   training: [nlZ dnlZ        ] = gp(hyp, inf, mean, cov, lik, x, y);
% prediction: [ymu ys2 nlZ post] = gp(hyp, inf, mean, cov, lik, x, y, xs);
%
% where:
%
%   hyp      column vector of hyperparameters
%   inf      function specifying the inference method 
%   cov      prior covariance function (see below)
%   mean     prior mean function
%   lik      likelihood function
%   x        n by D matrix of training inputs
%   y        column vector of length n of training targets
%   xs       ns by D matrix of test inputs
%   ys       column vector of length nn of test targets
%
%   nlZ      returned value of the negative log marginal likelihood
%   dnlZ     column vector of partial derivatives of the negative
%               log marginal likelihood w.r.t. each hyperparameter
%   ymu      column vector (of length ns) of predictive output means
%   ys2      quadratic matrix (of length ns) of predictive output covariances
%
%   post     struct representation of the (approximate) posterior
%            3rd output in training mode and 6th output in prediction mode
% 
% See also covFunctions.m, infMethods.m, likFunctions.m, meanFunctions.m.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2010-12-22
% KRZ- Added a global variable to store the test time of the algorithm.
if nargin<7 || nargin>8
  disp('Usage: [nlZ dnlZ          ] = gp(hyp, inf, mean, cov, lik, x, y);')
  disp('   or: [ymu ys2 nlZ post  ] = gp(hyp, inf, mean, cov, lik, x, y, xs);')
  return
end

if isempty(inf),  inf = @infExact; else                        % set default inf
  if iscell(inf), inf = inf{1}; end                      % cell input is allowed
  if ischar(inf), inf = str2func(inf); end        % convert into function handle
end
if isempty(mean), mean = {@meanZero}; end                     % set default mean
if ischar(mean) || isa(mean, 'function_handle'), mean = {mean}; end  % make cell
if isempty(cov), error('Covariance function cannot be empty'); end  % no default
if ischar(cov)  || isa(cov,  'function_handle'), cov  = {cov};  end  % make cell
cov1 = cov{1}; if isa(cov1, 'function_handle'), cov1 = func2str(cov1); end
if strcmp(cov1,'covFITC'); inf = @infFITC; end       % only one possible inf alg
if ~isempty(lik), error('Only [] supported for lik!'); end    
lik = @likGauss; 
D = size(x, 2); %will be used in the following fevals
if ~isfield(hyp,'mean'), hyp.mean = []; end        % check the hyp specification
if eval(feval(mean{:})) ~= numel(hyp.mean)
  error('Number of mean function hyperparameters disagree with mean function')
end
if ~isfield(hyp,'cov'), hyp.cov = []; end
  if eval(feval(cov{:})) ~= numel(hyp.cov)
  error('Number of cov function hyperparameters disagree with cov function')
end
if ~isfield(hyp,'lik'), hyp.lik = []; end
if eval(feval(lik)) ~= numel(hyp.lik)
  error('Number of lik function hyperparameters disagree with lik function')
end

if nargin==7                                     % if no test cases are provided
  [post, nlZ, dnlZ] = inf(hyp, mean, cov, lik, x, y);
  varargout = {nlZ, dnlZ, post};    % report -log marg lik, derivatives and post
else
  [post, nlZ] = inf(hyp, mean, cov, lik, x, y);
  alpha = post.alpha; L = post.L; sW = post.sW;
  Ltril = all(all(tril(L,-1)==0)) && ~all(all(L == 0));            % is L an upper triangular matrix?
  %ymu = zeros(ns,1); ys2 = zeros(ns, ns);   % allocate mem
  Ks  = feval(cov{:}, hyp.cov, x, xs);         % cross-covariances
  ms = feval(mean{:}, hyp.mean, xs);
  ymu = ms + Ks'*alpha;                       % predictive means
  if Ltril           % L is triangular => use Cholesky parameters (alpha,sW,L)
    Kss = feval(cov{:}, hyp.cov, xs);              % self-variance
    V  = L'\(repmat(sW,1,size(xs, 1)).*Ks);
    ys2 = Kss - V'*V;                       % predictive variances
  else                % L is not triangular => use alternative parametrisation
    Kss = feval(cov{:}, hyp.cov, xs, 'prior');              % self-variance
    ys2 = Kss + Ks'*L*Ks;                 % predictive variances
  end
  ys2 = max(ys2,0);   % remove numerical noise i.e. negative variances
  ys2 = ys2 + exp(2*hyp.lik)*eye(size(ys2, 1));
  varargout = {ymu, ys2, nlZ, post};        % assign output arguments
end
