function [cosine] = vectorCosine(x, y)
% VECTOR COSINE computes the cosine of the angle between the vectors x & y
% outputing a value between 1 and -1.
% The cosine is useful as a measure of similarity:
% 0 means the vector are orthogonal, or completely dissimilar
% +1 means the vectors are identical
% -1 means the vectors are exact opposites

[row, col] = size(x);

% if not a vector of singles, convert to a vector of singles
if(row == 1 || col == 1)
    a = single(x);
    b = single(y);
else
    a = single(reshape(x,[1,row*col]));
    b = single(reshape(y,[1,row*col]));
end

lengthX = dot(a,a);
lengthY = dot(b,b);

productXY = dot(a,b);

if lengthX == 0
    error('Error in vector cosine: first parameter has magnitude of zero.');
elseif lengthY == 0
    error('Error in vector cosine: second parameter has magnitude of zero.');
else
    cosine = productXY / (sqrt(lengthX) * sqrt(lengthY)); 
end

end