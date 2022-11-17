function h=fun_tiqu(x,L)
%提取出信息比特
%x为输入信号，L为信号的码元长度
for i=1:length(x)/L
    C=x((i-1)*L+1:i*L);
    t=sum(C);
    if t>=0
        h(i)=1;
    else
        h(i)=0;
    end          
end