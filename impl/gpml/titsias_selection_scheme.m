function [ cand, ind ] = titsias_selection_scheme( Candidates, Selected, k, X, y, sn2 )
%TITSIAS_SELECTION_SCHEME Summary of this function goes here
%   Detailed explanation goes here

M = size(Selected, 2);
n = size(X, 1); 
cand = 0;
ind = 0;
best = -Inf;
diagK = diag(k(X, X'));
for j = 1:size(Candidates, 2)
    U = [X(Candidates(j), :); X(Selected, :)];
    kUU = k(U, U');
    kU = k(U, X');
    snu2 = 1e-6;                              % hard coded inducing inputs noise
    Luu  = chol(kUU + snu2*eye(M+1));                       % Kuu + snu2*I = Luu'*Luu
    V  = Luu'\kU;                                     % V = inv(Luu')*Ku => V'*V = Q
    dg = diagK + sn2 - sum(V.*V,1)';      % D + sn2*eye(n) = diag(K) + sn2 - diag(Q)
    V  = V./repmat(sqrt(dg)',M+1,1);
    Lu = chol(eye(M+1) + V*V');
    r  = y./sqrt(dg);
    be = Lu'\(V*r);

    nlZ = sum(log(diag(Lu))) + (sum(log(dg)) + n*log(2*pi) + r'*r - be'*be)/2; 
    % Titsias bound
    iKuuKu = Luu'\kU;
    nlZ = nlZ + (sum(diagK) - trace(iKuuKu'*iKuuKu))/(2*sn2);
    
    if -nlZ > best
       best = -nlZ;
       cand = Candidates(j);
       ind = j;
    end
end
end

