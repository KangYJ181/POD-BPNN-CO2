% 1. This code is used to find the optimal q value.
% 2. Initial conditions:
%    (1) temperature: 22.5 degree celsius;
%    (2) pressure: 10 MPa;
%    (3) relative height of leak: 0.5.
% 3. Please load case_q.mat, FI_P.mat, U_P_ave.mat, NET_AP1.mat,..., NET_AP11.mat
% manually before running this code.

clc

% Initial conditions
T_initial=22.5;
P_initial=10;
Height_real=0.5;

% Length of the pressure response
M=10000;

% The pressure response
P=case_q(1:M,2);

% Initialization
height=[];  % Solution to the relative height of the leak
epsilon=[];  % Relative error
results=[];
num=0;  % Number of for cycles

% Change the q value and solve relative height of the leak.
% Take some divisors of M=10000 as values of q.
for q=[5,10,25,50,100,200,1250,10000]
    n=round(M/q);
    U_P_ave_q=U_P_ave(1:n:M);
    FI_P_q=FI_P(1:n:M,:);
    P_q=P(1:n:M);
    height=inverse_algorithm_q(T_initial,P_initial,P_q,...
                               U_P_ave_q,FI_P_q,NET_AP1,...
                               NET_AP2,NET_AP3,NET_AP4,...
                               NET_AP5,NET_AP6,NET_AP7,...
                               NET_AP8,NET_AP9,NET_AP10,...
                               NET_AP11);
    num=num+1;
    epsilon=abs(Height_real-height)/Height_real*100;
    results=[results;[height,epsilon]];
    [num,epsilon]
end