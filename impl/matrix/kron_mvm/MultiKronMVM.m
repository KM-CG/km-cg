function v = MultiKronMVM(A, Np, id, u)
    %if numel(size(A)) > 2 || size(A, 2) > 1, error('Please provide the Kronecker product as cell array.'); end
    D = size(A, 1); % How many Kronecker products?

    elements_missing = numel(id) < Np;
    
    if elements_missing
        v = zeros([Np, 1]); v(id) = u(:);
    else
        v = u;
    end
    
    for d = 1 : D
        B = A{d};
        v = (B * reshape(v, [size(B, 2) Np / size(B, 2)]))';
        %v = MatTimesTensor(B, v, d);
    end
    
    if elements_missing
        v = reshape(v(id), size(u));
    else
        v = reshape(v, size(u));
    end
return
