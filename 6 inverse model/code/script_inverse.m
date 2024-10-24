% 1.This code is used to solve the inverse problem.
% 2. Please load case_total.mat, FI_P.mat, U_P_ave.mat, NET_AP1.mat,..., NET_AP11.mat,
%InitialConditions.mat  manually before running this code.


clc
SNR=40;  % Signal to noise ratio[dB]
% Save results to the matrix.
NoisedP_and_Noise=[];
Height_and_error=[];

for k1=1
    % A pressure response of length q=100 is obtained.
    P=case_total(:,k1);
    P_q=P(1:100:10000);
    U_P_ave_q=U_P_ave(1:100:10000);
    FI_P_q=FI_P(1:100:10000,:);
    T_initial=InitialConditions(k1,1);
    P_initial=InitialConditions(k1,2);
    Height_real=InitialConditions(k1,3);

    % Add noises on pressure response ten times.
    for k2=1:1
    [NoisedP,Noise]=my_GaussianNoise(P_q,SNR);

    % Solve the relative height of leak.
    Height=inverse_algorithm_q(T_initial,P_initial,NoisedP,...
                               U_P_ave_q,FI_P_q,NET_AP1,...
                               NET_AP2,NET_AP3,NET_AP4,...
                               NET_AP5,NET_AP6,NET_AP7,...
                               NET_AP8,NET_AP9,NET_AP10,...
                               NET_AP11);
    if (Height>=0.0023)&&(Height<=0.9977)
        Height=Height;
    elseif Height<0.0023
        Height=0.0023;
    else
        Height=0.9977;
    end

    % Save the results.
    error=abs(Height-Height_real)/Height_real*100;
    Height_and_error=[Height_and_error,[Height;error]]
    NoisedP_and_Noise=[NoisedP_and_Noise,[NoisedP,Noise]];
    end
end
ave_error=mean(Height_and_error(2,:))