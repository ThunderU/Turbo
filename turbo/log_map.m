function L_all = log_map(rec_s,g,L_a,ind_dec)
% Log_MAP algorithm using straightforward method to compute branch metrics
% no approximation is used.
% Can be simplified to Max-Log-MAP by using approximation ln(e^x+e^y) = max(x,y).
% Input: rec_s: scaled received bits.
% rec_s = 0.5 * L_c * yk = ( 2 * a * rate * Eb/N0 ) * yk
% g: code generator for the component RSC code, in binary matrix form.
% L_a: a priori info. for the current decoder,
% scrambled version of extrinsic Inftyo. of the previous decoder.
% ind_dec: index of decoder. Either 1 or 2.
% Encoder 1 is assumed to be terminated, while encoder 2 is open.
% Output: L_all: log-likelihood ratio of the symbols. Complete information.
% Total number of bits: Inftyo. + tail
L_total = length(rec_s)/2;
[n,K] = size(g);
m = K - 1;
nstates = 2^m; % number of states in the trellis
% Set up the trellis
[next_out, next_state, last_out, last_state] = trellis(g);
Infty = 1e10;
% Initialization of Alpha
Alpha(1,1) = 0;
Alpha(1,2:nstates) = -Infty*ones(1,nstates-1); % first row of matrix Alpha
% Initialization of Beta
if ind_dec==1
    Beta(L_total,1) = 0;
    Beta(L_total,2:nstates) = -Infty*ones(1,nstates-1); % the last row of matrix Beta --yzh
elseif ind_dec==2
    Beta(L_total,1:nstates) = zeros(1,nstates); % the last row of matrix Beta --yzh
else
    fprintf('ind_dec is limited to 1 and 2!\n');
end
% Trace forward, compute Alpha
for k = 2:L_total+1
    for state2 = 1:nstates
        gamma = -Infty*ones(1,nstates);
        gamma(last_state(state2,1)) = (-rec_s(2*k-3)+rec_s(2*k-2)*last_out(state2,2))....
            -log(1+exp(L_a(k-1)));
        gamma(last_state(state2,2)) = (rec_s(2*k-3)+rec_s(2*k-2)*last_out(state2,4))....
            +L_a(k-1)-log(1+exp(L_a(k-1)));
        if(sum(exp(gamma+Alpha(k-1,:)))<1e-300)
            Alpha(k,state2)=-Infty;
        else
            Alpha(k,state2) = log( sum( exp( gamma+Alpha(k-1,:) ) ) );
        end
    end
    tempmax(k) = max(Alpha(k,:));
    Alpha(k,:) = Alpha(k,:) - tempmax(k);
end
% Trace backward, compute Beta
for k = L_total-1:-1:1
    for state1 = 1:nstates
        gamma = -Infty*ones(1,nstates);
        gamma(next_state(state1,1)) = (-rec_s(2*k+1)+rec_s(2*k+2)*next_out(state1,2))....
            -log(1+exp(L_a(k+1)));
        gamma(next_state(state1,2)) = (rec_s(2*k+1)+rec_s(2*k+2)*next_out(state1,4))....
            +L_a(k+1)-log(1+exp(L_a(k+1)));
        if(sum(exp(gamma+Beta(k+1,:)))<1e-300)
            Beta(k,state1)=-Infty;
        else
            Beta(k,state1) = log(sum(exp(gamma+Beta(k+1,:))));
        end
    end
    Beta(k,:) = Beta(k,:) - tempmax(k+1);
end

% Compute the soft output, log-likelihood ratio of symbols in the frame
for k = 1:L_total
    for state2 = 1:nstates
        gamma0 = (-rec_s(2*k-1)+rec_s(2*k)*last_out(state2,2))....
            -log(1+exp(L_a(k)));
        gamma1 = (rec_s(2*k-1)+rec_s(2*k)*last_out(state2,4))...
            +L_a(k)-log(1+exp(L_a(k)));
        temp0(state2) = exp(gamma0 + Alpha(k,last_state(state2,1)) + Beta(k,state2));
        temp1(state2) = exp(gamma1 + Alpha(k,last_state(state2,2)) + Beta(k,state2));
    end
    L_all(k) = log(sum(temp1)) - log(sum(temp0));
end