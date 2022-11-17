function en_output = encoderm( x, g, alpha, puncture )
% uses interleaver map 'alpha'
% if puncture = 1, unpunctured, produces a rate 1/3 output of fixed length
% if puncture = 0, punctured, produces a rate 1/2 output
% multiplexer chooses odd check bits from RSC1 and even check bits from RSC2
% determine the constraint length (K), memory (m)and number of information bits plus tail bits.
[n,K] = size(g);
m = K - 1;
L_info = length(x);
L_total = L_info + m;
% generate the codeword corresponding to the 1st RSC coder
% end = 1, perfectly terminated;
input = x;
output1 = rsc_encode(g,input,1);
% make a matrix with first row corresponing to info sequence
% second row corresponsing to RSC #1's check bits.
% third row corresponsing to RSC #2's check bits.
% y: unpuncture output of encoder; y(1,:) has m bits more than input bits
y(1,:) = output1(1:2:2*L_total);
y(2,:) = output1(2:2:2*L_total);
% interleave input to second encoder
for i = 1:L_total
    input1(1,i) = y(1,alpha(i)); %alpha--index of interleaver
end
%input has been interleaved. L_total bits already.end =-1?unterminated
output2 = rsc_encode(g, input1(1,1:L_total), -1 );
y(3,:) = output2(2:2:2*L_total);
% paralell to serial multiplex to get output vector
% puncture = 0: rate increase from 1/3 to 1/2;
% puncture = 1?unpunctured, rate = 1/3;
if puncture > 0 % unpunctured
    for i = 1:L_total
        for j = 1:3
            en_output(1,3*(i-1)+j) = y(j,i); % put the 3 bits of the same colomn to a sequential outputs
        end
    end
else % punctured into rate 1/2
    for i=1:L_total
        en_output(1,n*(i-1)+1) = y(1,i);
        if rem(i,2) % output check bits by turns
            en_output(1,n*i) = y(2,i); % odd check bits from RSC1
        else
            en_output(1,n*i) = y(3,i); % even check bits from RSC2
        end
    end
end
% antipodal modulation: +1/-1
en_output = 2 * en_output - ones(size(en_output));