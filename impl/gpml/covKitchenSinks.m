function K = covKitchenSinks( Omega, hyp, x, z, i )
%COVKITCHENSINKS Naive implementation of Random Kitchen sinks.
% Omega - M x D matrix of samples.

if nargin<3, K = '(D+1)'; return; end              % report number of parameters
if nargin<4, z = []; end                                   % make sure, z exists
xeqz = numel(z)==0; dg = strcmp(z,'diag') && numel(z)>0;        % determine mode

[n,D] = size(x);
ell = exp(hyp(1:D));                               % characteristic length scale
sf2 = exp(2*hyp(D+1));                                         % signal variance
x   = x./repmat(ell', [n, 1]);
M   = Omega * x';
Phi = [cos(M); sin(M)];

if dg                                                               % vector kxx  
  K   = sf2 * (Phi' * Phi) / size(Omega, 1);
  K   = diag(K);
else
  if xeqz                                                 % symmetric matrix Kxx
    K   = sf2 * (Phi' * Phi) / size(Omega, 1);
  else                                                   % cross covariances Kxz
    z   = z./repmat(ell', [size(z, 1), 1]);
    M   = Omega * z';
    Phiz= [cos(M); sin(M)];
    K = sf2 * (Phi' * Phiz) / size(Omega, 1);
  end
end

if nargin>4                                                        % derivatives
    error('Derivatives not supported.')
end
end

