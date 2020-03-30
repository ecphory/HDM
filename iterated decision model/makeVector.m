function [ vector ] = makeVector(n)
%makeVector
% Generate a vector of n random values
% selected from an normal distribution
% with a mean of 0 and a variance of 1/n
% (i.e. a standard deviation of sqrt(1/n))

sd   = sqrt(1/n);
mean = 0;

vector = mean + sd*randn(n,1);

end