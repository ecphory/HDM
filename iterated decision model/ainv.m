function [ vrotce ] = ainv( vector )
%AINV returns the approximate inverse of a vector
%   vrotce is vector with all elements reversed save the first

n = length(vector);

vrotce = vector;

for index = 2:n
    vrotce(index) = vector(n + 2 - index);
end

end