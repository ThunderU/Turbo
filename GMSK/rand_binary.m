function [data_binary,data_binary1]=rand_binary(data_len);
%�������һ��������������Ϊ�����õ���Ϣ����
%**************************************************************************
%data           ���г���
%data_binary    �����Ķ���������
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