function phase = gmsk_phase(data,data_len,sample_number,qt)
%计算调制相位
%**************************************************************************
%data           输入信号序列
%data_len       序列长度
%sample_nubmer  采样个数
%qt             q(t)
%**************************************************************************

%--------------------------------------------------------------------------
%累积相位
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
%即时相位
defai = zeros(1,data_len*sample_number);
defai(1:sample_number) = pi*(data(2) * qt(1:sample_number) + data(1) * qt(sample_number+1:2*sample_number));
for i = 2:data_len-1
    defai((i-1)*sample_number+1:i*sample_number) = pi*(data(i-1) * qt(2*sample_number+1:3*sample_number) + data(i) * qt(sample_number+1:2*sample_number) + data(i+1) * qt(1:sample_number));
end
defai((data_len-1)*sample_number+1:data_len*sample_number) = pi*(data(data_len-1) * qt(2*sample_number+1:3*sample_number) + data(data_len) * qt(sample_number+1:2*sample_number));
%**************************************************************************

%--------------------------------------------------------------------------
%调制相位
phase = defai + theta1;
%**************************************************************************