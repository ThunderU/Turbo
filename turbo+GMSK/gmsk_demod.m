function demod_data = gmsk_demod(signal_mod1)
    sample_number = 8;          %采样个数
    Rb = 24000;                 %码元速率
    fc = 96000;                 %载波频率
    alpha = 0.3;               %BbTb值
    multi = fc/Rb;
    data_len = length(signal_mod1)/multi/sample_number;
    Fs=fc*sample_number;

    cospart = real(signal_mod1);
    sinpart = imag(signal_mod1);
    
    delay = [zeros(1,2*multi*sample_number) cospart(1:length(signal_mod1)-2*multi*sample_number)];
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