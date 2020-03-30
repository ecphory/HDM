function [ echo ] = ccorr(probe,trace)
%CCORR computes the circular correlation of vectors probe and trace
%   by taking the approximate inverse (ainv) of the probe and
%   convolving it (cconv) with the trace.
%   If trace is the circular convolution (cconv)
%   of probe and another vector, then the output of this function 
%   will be a likeness of that other vector

n = length(trace);

pebor = ainv(probe);
echo = cconv(pebor,trace,n);

end