function [ a ] = realNorm( a )
%realNorm normalizes a vector to a mean of 0 and Euclidean length of 1

a = a - mean(a);
a = a / sqrt(dot(a,a));

end

