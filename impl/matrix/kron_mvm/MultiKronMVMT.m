function v = MultiKronMVMT(A, Np, id, u)
    D = size(A, 1); % How many Kronecker products?
    
    if numel(id) < Np
        v = zeros([Np, 1]); v(id) = u(:);
    else
        v = u;
        end
    
    for d = 1 : D
        B = A{d};
        v = (B' * reshape(v, [size(B, 2) Np / size(B, 2)]))';
    end
    v = reshape(v(id), size(u));
return
