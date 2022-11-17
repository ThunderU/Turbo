%GMSK调制，二比特差分解调

function signal_mod1 = gmsk_mod(data)
%--------------------------------------------------------------------------
%参数设置
data_len = length(data);          %码元个数
sample_number = 8;          %采样个数
Rb = 24000;                 %码元速率
fc = 96000;                 %载波频率
alpha = 0.3;               %BbTb值
%**************************************************************************

%--------------------------------------------------------------------------
%随机产生传输信号
% data = rand_binary(data_len);
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

% signal_mod_i=signal_mod_i;
% signal_mod_q=signal_mod_q;
signal_mod1=signal_mod_i+1i*signal_mod_q;
%--------------------------------------------------------------------------