  %GMSK调制，二比特差分解调
clear all
close all

%--------------------------------------------------------------------------
%参数设置
data_len = 200000;          %码元个数
sample_number = 8;          %采样个数
Rb = 24000;                 %码元速率
fc = 96000;                 %载波频率
alpha = 0.25;               %BbTb值
%**************************************************************************

%--------------------------------------------------------------------------
%随机产生传输信号
data = rand_binary(data_len);
%**************************************************************************

%--------------------------------------------------------------------------
%差分编码
data_diff = difference(data);
%**************************************************************************

%--------------------------------------------------------------------------
%参数设置
[signal_out,I_out,Q_out,phase,gt,qt] = mod_gmsk(data_diff,data_len,sample_number,Rb,alpha);
%**************************************************************************

%--------------------------------------------------------------------------
%中频搬移
multi = fc/Rb;
I_temp=interp(I_out,multi);
Q_temp=interp(Q_out,multi);

Fs=fc*sample_number;
t=1/Fs:1/Fs:length(I_temp)*1/Fs;
signal_i=I_temp.*cos(2*pi*fc*t);
signal_q=Q_temp.*sin(2*pi*fc*t);
signal_mod_i=I_temp.*cos(2*pi*fc*t)-Q_temp.*sin(2*pi*fc*t);
signal_mod_q=I_temp.*sin(2*pi*fc*t)+Q_temp.*cos(2*pi*fc*t);
%**************************************************************************

signal_mod_i=signal_mod_i;
signal_mod_q=signal_mod_q;
signal_mod1=signal_mod_i+j*signal_mod_q;

[p1,f1]=periodogram(signal_mod1);
plot(f1,10*log10(p1));
xlabel('归一化频率')
ylabel('功率谱密度/dB')
title('GMSK功率谱密度')


%--------------------------------------------------------------------------
%加噪声
% for SNR = 0:16
    
% signal_mod1 = awgn(signal_mod_i+j*signal_mod_q,SNR);

% noise=randn(1,length(signal_mod_i));
% b=5/50;c=6/50;d=19/50;e=20/50;
% F1=[0,b,c,d,e,1];
% A1=[0,0,1,1,0,0];
% bpf=firls(300,F1,A1);
% 
% noise1=conv(noise,bpf);
% noise_out=noise1(150+1:150+length(signal_mod_i));
% 
% Power_noise=sum(noise_out.^2);
% Power_signal_i=sum(signal_mod_i.^2);
% Power_signal_q=sum(signal_mod_q.^2);
% Coe_SNR_i=sqrt(1*Power_signal_i/Power_noise/(10^(SNR/10)));
% Coe_SNR_q=sqrt(1*Power_signal_q/Power_noise/(10^(SNR/10)));
% 
% signal_mod1 = (signal_mod_i+(Coe_SNR_i)*noise_out) + j * (signal_mod_q+(Coe_SNR_q)*noise_out);    %加入噪声以后的信号

    %--------------------------------------------------------------------------

    EbN0 = 0:0.5:16;
    beri = zeros(1,length(EbN0));
 for i = 1:length(EbN0)
     for k=1:100
        signal_mod1 = awgn(signal_mod_i+1i*signal_mod_q,EbN0(i));
    
        %--------------------------------------------------------------------------
        %二比特差分
        cospart = real(signal_mod1);
        sinpart = imag(signal_mod1);
        
        delay = [zeros(1,2*multi*sample_number) cospart(1:data_len*sample_number*multi-2*multi*sample_number)];
        cha = -delay .* cospart;
        
        N=300;                                              % 滤波器的阶数为(N+1)  
        F=[0,fc-1000,fc+1000,Fs/2]*2/Fs;
        A=[1,1,0,0];
        lpf=firls(N,F,A);
        
        dem1 = conv(cha,lpf);
        dem = dem1(N/2+1:N/2+length(cha));
        %**************************************************************************
        
        %--------------------------------------------------------------------------
        %抽样判决
        for j = 2:data_len
            demod_data(j) = dem(j*multi*sample_number)>-0.1;
        end
        demod_data = 2*demod_data-1;
        %**************************************************************************
        
        %--------------------------------------------------------------------------
        %计算误码率
        [num,ber(i)]=symerr(demod_data(3:data_len),data(3:data_len));
        beri(i) = beri(i)+ber(i);
    end
    beri(i)=beri(i)/100;
    fprintf('%d',i); 
end
%**************************************************************************

%--------------------------------------------------------------------------
%误码率曲线
semilogy(EbN0,ber,'r*-');
title('GMSK误码率');
set(gca,'FontSize',24);
xlabel('SNR/dB','FontSize',16);
ylabel('BER','FontSize',16);
grid on;
%**************************************************************************