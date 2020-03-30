function [ chunk ] = hdmUOG( percepts, p, left)
%hdmUOG is a variant of hemUOG designed for Holographic Declarative Memory
%   Unlike hemUOG, hdmUOG is designed to ONLY include UOGs
%   that include the designated item p, where p is an index
%
%   In the HDM model, p is the placeholder vector.
%
%   UOG or Unconstrained Open n-Grams are all ordered sequences of
%   characters/vectors in a string, up to grams of numOfItems characters,
%   and including sequences that "skip" characters, i.e.,
%   sequences are not restricted by adjacency.
%
%   e.g., all UOG of the sequence 'abc' will be:
%   a, b, c, ab, ac, bc, abc.
%
%   percepts: matrix of environmental vectors of dimensions [numOfItems,N]
%             where: numOfItems is the number of vectors
%                    N is the dimensionality of each vector
%
%   chunk: sum of all UOG of the environmental vectors in percepts
%
%   p: the index of the placeholder vector, i.e., the item that must be
%      included in all skip-grams / UOGs
%          
%   left: permutation for indicating that a vector is the left operand of
%         a non-commutative circular convolution. By default, the function
%         uses no permutation (i.e., the numbers 1 to N in ascending order)

[numOfItems,N] = size(percepts);

if nargin < 3
    left = 1:N;
end

chunk = zeros(1,N);
sum   = zeros(1,N);
for i=1:numOfItems
    if i == 1
        sum = percepts(i,:); 
    elseif (i > 1) && (i < p)
        leftOperand = chunk + sum;
        chunk = chunk + cconv(leftOperand(left),percepts(i,:),N);
        sum = sum + percepts(i,:);
    elseif i == p % force all skip grams to include item p
        leftOperand = chunk + sum;
        chunk = cconv(leftOperand(left),percepts(i,:),N);
        sum = percepts(i,:);
    else % i > p, i > 1
        leftOperand = chunk + sum;
        chunk = chunk + cconv(leftOperand(left),percepts(i,:),N);
    end
end

end