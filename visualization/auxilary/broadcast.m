function c = broadcast(f,a,b,def)
% broadcasts a function over different length using 'def' as default value
    if nargin < 4, def = 0; end
    num_steps = length(b);
    if num_steps > length(a)
        c = [a, def * ones(1, num_steps - length(a))];
    else
        c = a;
    end
    c(1:num_steps) = f(c(1:num_steps), b);
end

