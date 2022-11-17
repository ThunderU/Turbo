function [next_out, next_state, last_out, last_state] = trellis(g)
% set up the trellis given code generator g
% g given in binary matrix form. e.g. g = [ 1 1 1; 1 0 1 ];
% next_out(i,1:2): trellis next_out (systematic bit; parity bit) when input = 0, state = i; next_out(i,j)
%= -1 or 1
% next_out(i,3:4): trellis next_out (systematic bit; parity bit) when input = 1, state = i;
% next_state(i,1): next state when input = 0, state = i; next_state(i,i) = 1,...2^m
% next_state(i,2): next state when input = 1, state = i;
% last_out(i,1:2): trellis last_out (systematic bit; parity bit) when input = 0, state = i; last_out(i,j)
%= -1 or 1
% last_out(i,3:4): trellis last_out (systematic bit; parity bit) when input = 1, state = i;
% last_state(i,1): previous state that comes to state i when info. bit = 0;
% last_state(i,2): previous state that comes to state i when info. bit = 1;
[n,K] = size(g);
m = K - 1;
max_state = 2^m;
% set up next_out and next_state matrices for systematic code
for state=1:max_state
    state_vector = bin_state( state-1, m ); % matrix state_vector is of max_state rows and m columns
    % when receive a 0
    d_k = 0;
    a_k = rem( g(1,:)*[0 state_vector]', 2 );
    [out_0, state_0] = encode_bit(g, a_k, state_vector);
    out_0(1) = 0;
    % when receive a 1
    d_k = 1;
    a_k = rem( g(1,:)*[1 state_vector]', 2 );
    [out_1, state_1] = encode_bit(g, a_k, state_vector);
    out_1(1) = 1;
    next_out(state,:) = [out_0 out_1]; % BPSK? Each row has two possible outputs(according to input 1
    %or 0) --yzh
    next_state(state,:) = [(int_state(state_0)+1) (int_state(state_1)+1)]; % 2 next state for current state
    %according to input --yzh
end

% find out which two previous states can come to present state
last_state = zeros(max_state,2);
for bit=0:1
    for state=1:max_state
        last_state(next_state(state,bit+1), bit+1)=state; % row number is the next_state, column is the input bit --yzh
        last_out(next_state(state, bit+1), bit*2+1:bit*2+2) ... % row is the next_state value --yzh
            = next_out(state, bit*2+1:bit*2+2); % next_out is the output of current state with input 0 or 1 -yzh
    end
end
