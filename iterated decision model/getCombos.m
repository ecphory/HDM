function [ combo ] = getCombos( percepts, placeholder, cLambda, p, left)
%getCombos
%   is a function for getting all ordered sequences of vectors in percepts
%   and summing them together to create a vector referred to here as
%   'combos'. This is used for probing memory in BEAGLE and DSHM.
%
%   percepts: matrix of environmental vectors of dimensions [numOfItems,N]
%             where: numOfItems is the number of vectors
%                    N is the dimensionality of each vector
%
%   combo: sum of all sequences of the environmental vectors in
%          percepts, up to sequences, "chunks", of cLambda vectors
%          with the p-th vector replaced with the placeholder
%          (so, uh, it doesn't matter what you put there)
%
%   placeholder: placeholder vector, the "key" to all memory
%
%   cLambda: maximum chunk size
%
%   p: index of the placeholder
%
%   left: permutation for indicating that a vector is the left operand of
%         a non-commutative circular convolution

[numOfItems,N] = size(percepts);

if nargin < 5
    left = 1:N;
end

combo = zeros(1,N);

for j=1:p             % iterate on starting point of chunk
    for i=j:(min(j + cLambda - 1,numOfItems)) % iterate within chunk
        if i == p
            thisPercept = placeholder;
        else
            thisPercept = percepts(i,:);
        end
        
        if i == j
            chunk = thisPercept;
        elseif i < p
            chunk = cconv(chunk(left),thisPercept,N);
        else
            chunk = cconv(chunk(left),thisPercept,N);
            combo = combo + chunk;
        end
    end
end


end