function bin_state = bin_state( int_state, m )
% converts an vector of integer into a matrix; the i-th row is the binary form
% of m bits for the i-th integer
for j = 1:length( int_state ) % length(int_state)?=max_state? --yzh
    for i = m:-1:1
        state(j,m-i+1) = fix( int_state(j)/ (2^(i-1)) ); % fix(X) rounds the elements of X to the nearest
        %integers towards zero.
        int_state(j) = int_state(j) - state(j,m-i+1)*2^(i-1); % remain of mod 2^(i-1), the leftmost bit
        %is most significant -yzh
    end
end
bin_state = state;