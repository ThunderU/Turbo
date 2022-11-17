function y = rsc_encode(g, x, terminated)
% encodes a block of data x (0/1)with a recursive systematic convolutional code
% with generator vectors in g, and returns the output in y (0/1).
% if terminated>0, the trellis is perfectly terminated
% if terminated<0, it is left unterminated;
% determine the constraint length (K), memory (m), and rate (1/n)and number of information bits.
[n,K] = size(g);
m = K - 1;
if terminated>0
    L_info = length(x); % L_info: lenght of information sequence
    L_total = L_info + m; % L_total:m additional bits is used to terminate
else
    L_total = length(x);
    L_info = L_total - m; % see the sequence for untermated in function encoderm for reason. length of x is
    L_total;
end
state = zeros(1,m); % initialize the state vector

% generate the codeword
for i = 1:L_total
    if terminated<0 || (terminated>0 && i<=L_info)
        d_k = x(1,i); % d_k: information sequence
    elseif terminated>0 && i>L_info % terminate the trellis
        d_k = rem( g(1,2:K)*state', 2 ); %??a
    end
    a_k = rem( g(1,:)*[d_k state]', 2 ); % a_k: the bit to be put into the register
    [output_bits, state] = encode_bit(g, a_k, state);
    output_bits(1,1) = d_k; % since systematic, first output is input bit
    y(n*(i-1)+1:n*i) = output_bits; % n output bits for 1 input bit(recursiv encoder)
end
