function phase = gmsk_phase(data,data_len,sample_number,qt)
%���������λ
%**************************************************************************
%data           �����ź�����
%data_len       ���г���
%sample_nubmer  ��������
%qt             q(t)
%**************************************************************************

%--------------------------------------------------------------------------
%�ۻ���λ
theta = zeros(1,data_len);
for i = 3:data_len
    theta(i) = theta(i-1) + pi/2 * data(i-2);
end
theta1 = zeros(1,data_len*sample_number);
for i = 1:sample_number
    theta1(i:sample_number:data_len*sample_number) = theta;
end
%**************************************************************************

%--------------------------------------------------------------------------
%��ʱ��λ
defai = zeros(1,data_len*sample_number);
defai(1:sample_number) = pi*(data(2) * qt(1:sample_number) + data(1) * qt(sample_number+1:2*sample_number));
for i = 2:data_len-1
    defai((i-1)*sample_number+1:i*sample_number) = pi*(data(i-1) * qt(2*sample_number+1:3*sample_number) + data(i) * qt(sample_number+1:2*sample_number) + data(i+1) * qt(1:sample_number));
end
defai((data_len-1)*sample_number+1:data_len*sample_number) = pi*(data(data_len-1) * qt(2*sample_number+1:3*sample_number) + data(data_len) * qt(sample_number+1:2*sample_number));
%**************************************************************************

%--------------------------------------------------------------------------
%������λ
phase = defai + theta1;
%**************************************************************************