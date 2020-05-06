function [ m, mi, ma] = compute_statistics( v, m, mi, ma, v1, v2 )
% computes sequential mean, minimum and maximum of arrays of different length
    delta = broadcast(@(a, b) b - a, m ./ v1, v ./ v2);
    m  = broadcast(@(a, b) a + b, m, delta(1:length(v2)));
    mi = broadcast(@min, mi, v, Inf);
    ma = broadcast(@max, ma, v);
end