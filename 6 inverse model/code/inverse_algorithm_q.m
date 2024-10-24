% 功能：反演破口高度
% T0：初始温度
% P0：初始压力
% CASE_q：压力响应
% U_P_ave_q：压力平均向量
% FI_P_q：降阶基底
% NET_AP：训练好的降阶系数神经网络

function best_height=inverse_algorithm_q(T0,P0,CASE_q,U_P_ave_q,FI_P_q,...
                                         NET_AP1,NET_AP2,NET_AP3,...
                                         NET_AP4,NET_AP5,NET_AP6,...
                                         NET_AP7,NET_AP8,NET_AP9,...
                                         NET_AP10,NET_AP11)
CASE_wavy=CASE_q-U_P_ave_q;
error1=[]; error2=[]; error3=[]; error4=[]; % 初始化反演误差
H_guess_1=0.1:0.1:0.9;

    % 在0.1~0.9之间以0.1为间隔筛选出最接近的高度
    for k1=1:9
        A=[NET_AP1([T0;P0;H_guess_1(k1)]);NET_AP2([T0;P0;H_guess_1(k1)]);...
           NET_AP3([T0;P0;H_guess_1(k1)]);NET_AP4([T0;P0;H_guess_1(k1)]);...
           NET_AP5([T0;P0;H_guess_1(k1)]);NET_AP6([T0;P0;H_guess_1(k1)]);...
           NET_AP7([T0;P0;H_guess_1(k1)]);NET_AP8([T0;P0;H_guess_1(k1)]);...
           NET_AP9([T0;P0;H_guess_1(k1)]);NET_AP10([T0;P0;H_guess_1(k1)]);...
           NET_AP11([T0;P0;H_guess_1(k1)])];
        dCASE=CASE_wavy-FI_P_q*A;
        error1(k1)=sum(dCASE.*dCASE);
    end
    best_height_1=H_guess_1(find(error1==min(error1)));

    H_guess_2=(-0.1:0.01:0.1)+best_height_1;
    for k2=1:21
        A=[NET_AP1([T0;P0;H_guess_2(k2)]);NET_AP2([T0;P0;H_guess_2(k2)]);...
           NET_AP3([T0;P0;H_guess_2(k2)]);NET_AP4([T0;P0;H_guess_2(k2)]);...
           NET_AP5([T0;P0;H_guess_2(k2)]);NET_AP6([T0;P0;H_guess_2(k2)]);...
           NET_AP7([T0;P0;H_guess_2(k2)]);NET_AP8([T0;P0;H_guess_2(k2)]);...
           NET_AP9([T0;P0;H_guess_2(k2)]);NET_AP10([T0;P0;H_guess_2(k2)]);...
           NET_AP11([T0;P0;H_guess_2(k2)])];
        dCASE=CASE_wavy-FI_P_q*A;
        error2(k2)=sum(dCASE.*dCASE);
    end
    best_height_2=H_guess_2(find(error2==min(error2)));

    H_guess_3=(-0.01:0.001:0.01)+best_height_2;
    for k3=1:21
        A=[NET_AP1([T0;P0;H_guess_3(k3)]);NET_AP2([T0;P0;H_guess_3(k3)]);...
           NET_AP3([T0;P0;H_guess_3(k3)]);NET_AP4([T0;P0;H_guess_3(k3)]);...
           NET_AP5([T0;P0;H_guess_3(k3)]);NET_AP6([T0;P0;H_guess_3(k3)]);...
           NET_AP7([T0;P0;H_guess_3(k3)]);NET_AP8([T0;P0;H_guess_3(k3)]);...
           NET_AP9([T0;P0;H_guess_3(k3)]);NET_AP10([T0;P0;H_guess_3(k3)]);...
           NET_AP11([T0;P0;H_guess_3(k3)])];
        dCASE=CASE_wavy-FI_P_q*A;
        error3(k3)=sum(dCASE.*dCASE);
    end
    best_height_3=H_guess_3(find(error3==min(error3)));

    H_guess_4=(-0.001:0.0001:0.001)+best_height_3;
    for k4=1:21
        A=[NET_AP1([T0;P0;H_guess_4(k4)]);NET_AP2([T0;P0;H_guess_4(k4)]);...
           NET_AP3([T0;P0;H_guess_4(k4)]);NET_AP4([T0;P0;H_guess_4(k4)]);...
           NET_AP5([T0;P0;H_guess_4(k4)]);NET_AP6([T0;P0;H_guess_4(k4)]);...
           NET_AP7([T0;P0;H_guess_4(k4)]);NET_AP8([T0;P0;H_guess_4(k4)]);...
           NET_AP9([T0;P0;H_guess_4(k4)]);NET_AP10([T0;P0;H_guess_4(k4)]);...
           NET_AP11([T0;P0;H_guess_4(k4)])];
        dCASE=CASE_wavy-FI_P_q*A;
        error4(k4)=sum(dCASE.*dCASE);
    end
    best_height=H_guess_4(find(error4==min(error4)));

end