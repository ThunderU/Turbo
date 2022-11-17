close all;
clear all;
%%
%系统参数设计
T_start=0;%开始时间
T_stop=1000;%截止时间
T=T_stop-T_start;%仿真持续时间

rs=10;%传输速率
NumBits=20000;%传输bit数

fc=2.5*NumBits/T;%载波频率
L=100;%码元长度

N_sample=NumBits*L;% 采样点数
T_sample=T/N_sample;%采样间隔
f_sample=1/T_sample;%采样频率

f_res=f_sample/N_sample;%频率分辨率

Tb=T/NumBits;%码元周期，单个码元持续时间

%载波
n=1:N_sample;
c1=cos(2*pi*fc*n*T_sample);
c2=sin(2*pi*fc*n*T_sample);

%%
%MSK信号的调制

%随机产生原始01信息比特
a=rand(1,NumBits)>0.5;

%差分编码
a_cf=zeros(1,NumBits);
a_cf(1)=a(1);
for i=2:NumBits
    a_cf(i)=~xor(a(i),a_cf(i-1));%进行同或操作
end


%串并变换
A=reshape(a_cf,2,NumBits/2);
a_I=A(1,:);
a_Q=A(2,:);

%0转换为-1，1转换为1
a_I=a_I*2-1;
a_Q=a_Q*2-1;


%与余弦相乘进行脉冲成型,得到基带信号
n=1-L:L;
e_I=[];e_Q=[];
for i=1:NumBits/2
    e_I=[e_I,a_I(i)*cos((pi*n*T_sample)/(2*Tb))];
    e_Q=[e_Q,a_Q(i)*cos((pi*n*T_sample)/(2*Tb))];
end

% Q delay
number_delay=length(n)/2;
e_Q1=[zeros(1,number_delay),e_Q(1:length(e_Q)-number_delay)];
scatter(e_I,e_Q1)

%与载波相乘
sI=e_I.*c1;%I支路信号
sQ=e_Q.*c2;%Q支路信号
s=sI-sQ;%MSK调制信号

window=boxcar(length(s));                 %矩形窗 
nfft=1024; 
[Pxx,w]=periodogram(s);
figure;
plot(w,10*log(Pxx));
xlabel('归一化频率')
ylabel('功率谱密度/dB')
title('MSK功率谱密度')
ylim([-100 50]);
xlim([0 0.5]);

%%%信道%%%%%%

%设置信噪比
ebn0 = 0:8;
snr = ebn0 - 10*log10(0.5*16);
BER_test = zeros(1,length(snr));
err = zeros(1,length(BER_test));
%瑞利信道
for m = 1:length(snr)
    for k=1:100
        %线性高斯白噪声信道
        y_n = awgn(s,snr(m),'measured');
        
        %MSK信号的解调
    
        %通过带通滤波器，去除噪声
        [b1,a1]=butter(4,[0.02,0.06]);%计算带通滤波器的H(z)系数
        y=filtfilt(b1,a1,y_n);%对信号y_i进行滤波，得到信号y
    
        %与恢复载波相乘
        x1_I=y.*c1;
        x1_Q=-y.*c2;
    
        %通过低通滤波器，分离出两路信号
        [b2,a2]=butter(4,0.02,'low');%计算低通滤波器H(z)系数
        x2_I=filtfilt(b2,a2,x1_I);%对信号x1_I进行滤波，得到信号x2_I
        x2_Q=filtfilt(b2,a2,x1_Q);%对信号x1_Q进行滤波，得到信号x2_Q
    
        %提取出经过差分编码的01信息比特
        d_I=fun_tiqu(x2_I,2*L);
        d_Q=fun_tiqu(x2_Q,2*L);
    
        %并串变换,
        d0(1:2:NumBits)=d_I;
        d0(2:2:NumBits)=d_Q;
    
        %差分解码，恢复出原始01信息比特
        d=zeros(1,NumBits);
        d(1)=d0(1);
        for i=2:NumBits
            d(i)=~xor(d0(i-1),d0(i));%进行同或操作
        end
    
        %%
        %误码率计算   
        err(m)=err(m)+length(find(d~=a));%计算解调信号中错误码元个数
        
    end
    BER_test(m)=err(m)/NumBits/100;
end

%%
%图像
% figure(1);
% subplot(4,2,1:2);stem(a);title('原始01信息比特序列a');axis([1,20,0,1]);
% subplot(4,2,3:4);stem(a_cf);title('a经过差分编码后得到a_cf');axis([1,20,0,1]);
% subplot(4,2,5);stem(a_I);title('I支路-1,1比特序列a_I');axis([1,10,-1,1]);
% subplot(4,2,6);stem(a_Q);title('Q支路-1,1比特序列a_Q');axis([1,10,-1,1]);
% subplot(211);plot(e_I);title('I支路进行脉冲成型后得到e_I');axis([1,2000,-1.5,1.5]);
% subplot(212);plot(e_Q);title('Q支路进行脉冲成型后得到e_Q');axis([1,2000,-1.5,1.5]);
% 
% figure(2);
% subplot(2,2,1);plot(sI);title('I支路调制信号');axis([1,2000,-2,2]);
% subplot(2,2,2);plot(sQ);title('Q支路调制信号');axis([1,2000,-2,2]);
% subplot(2,2,3:4);plot(s);title('MSK调制信号');axis([1,4000,-2,2]);
% 
% figure(3);
% subplot(2,1,1);plot(s);title('MSK调制信号');axis([1,4000,-2,2]);
% subplot(2,1,2);plot(y_n);title('MSK调制信号s通过高斯信道后的信号y_n');axis([1,4000,-2,2]);
% 
% 
% figure(4);
% subplot(3,2,1:2);plot(y);title('y_i通过带通滤波器后的信号y');axis([0,4000,-2,2]);
% subplot(3,2,3);plot(x1_I);title('y与恢复载波c1相乘后的信号x1_I');axis([0,2000,-2,2]);
% subplot(3,2,4);plot(x1_Q);title('y与恢复载波c2相乘后的信号x1_Q');axis([0,2000,-2,2]);
% subplot(3,2,5);plot(x2_I);hold on;plot(e_I);
% legend('x1_I通过低通滤波器后的信号x2_I','I支路基带信号e_I');axis([0,2000,-2,2]);
% subplot(3,2,6);plot(x2_Q);hold on;plot(e_Q);
% legend('x1_Q通过低通滤波器后的信号x2_Q','Q支路基带信号e_Q');axis([0,2000,-2,2]);
% 
% figure(5);
% subplot(3,2,1);stem(d_I);
% title('x2_I经过判决后得到的01序列d_I');axis([1,10,0,1]);
% subplot(3,2,2);stem(d_Q);
% title('x2_Q经过判决后得到的01序列d_Q');axis([1,10,0,1]);
% subplot(3,2,3:4);stem(d0);hold on;stem(a_cf,'*');
% legend('d_I和d_Q并串变换后得到的01序列d0','原始经差分编码后的01序列a_cf');axis([1,20,0,1.5]);
% subplot(3,2,5:6);stem(d);hold on;stem(a,'*');
% legend('d0经差分解码后得到的解调01信息比特d','原始01信息比特a');axis([1,20,0,1.5]);

figure(6);
semilogy(ebn0,BER_test,'-o');
xlabel('比特信噪比');
ylabel('误码率');
title('MSK误码率曲线');
legend('实验曲线');
grid on