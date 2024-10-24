% 作用：产生额定信噪比的随机高斯噪声
% x_input：输入信号（列信号）
% snr：信噪比[dB]
% x_output：加噪信号（列信号）
function [x_output,noise]=my_GaussianNoise(x_input,snr)
L=length(x_input);                                                         % 信号长度
Ps=sum(x_input.^2)/L;                                                      % 信号功率
Pn=Ps/(10^(snr/10));                                                       % 噪声功率
noise=sqrt(Pn).*randn(L,1);                                                % 加噪信号
x_output=x_input+noise;
end