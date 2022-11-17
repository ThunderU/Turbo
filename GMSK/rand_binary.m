function [data_binary,data_binary1]=rand_binary(data_len);
%随机产生一个二进制序列作为仿真用的消息序列
%**************************************************************************
%data           序列长度
%data_binary    产生的二进制序列
%**************************************************************************

%--------------------------------------------------------------------------
data1=randn(1,data_len);
data_binary=zeros(1,data_len);
data_binary1=zeros(1,data_len);
for i=1:data_len
    if data1(i)<0
       data_binary(i)=-1;
       data_binary1(i)=0;
   else
       data_binary(i)=1;
       data_binary1(i)=1;
    end
end
%**************************************************************************