function [output, state] = encode_bit(g, input, state)
% This function takes as an input a single bit to be encoded, as well as
% the coeficients of the generator polynomials and the current state vector.
% It returns as output n encoded data bits, where 1/n is the code rate.
% the rate is 1/n ?k is the constraint length ?m is the amount of memory
[n,k] = size(g); 
m = k-1; %//i>j>a
% determine the next output bit
for i=1:n
output(i) = g(i,1)*input; % the first bit a_k's contribution to output
    for j = 2:k
        output(i) = xor(output(i),g(i,j)*state(j-1)); % a_(k-j)'s contribution to output
    end
end
state = [input, state(1:m-1)]; % shift one bit