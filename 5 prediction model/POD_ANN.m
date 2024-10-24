close
%输入待预测工况
T=29.1; P=10.31; H=0.4520;
%用神经网络还原系数
AP1=NET_AP1([T;P;H]);
AP2=NET_AP2([T;P;H]);
AP3=NET_AP3([T;P;H]);
AP4=NET_AP4([T;P;H]);
AP5=NET_AP5([T;P;H]);
AP6=NET_AP6([T;P;H]);
AP7=NET_AP7([T;P;H]);
AP8=NET_AP8([T;P;H]);
AP9=NET_AP9([T;P;H]);
AP10=NET_AP10([T;P;H]);
AP11=NET_AP11([T;P;H]);

AH1=NET_AH1([T;P;H]);
AH2=NET_AH2([T;P;H]);
AH3=NET_AH3([T;P;H]);
AH4=NET_AH4([T;P;H]);
AH5=NET_AH5([T;P;H]);
AH6=NET_AH6([T;P;H]);
AH7=NET_AH7([T;P;H]);
AH8=NET_AH8([T;P;H]);
AH9=NET_AH9([T;P;H]);
%还原对应工况下的压强/高度系数列向量
AP=[AP1;AP2;AP3;AP4;AP5;AP6;AP7;AP8;AP9;AP10;AP11];
AH=[AH1;AH2;AH3;AH4;AH5;AH6;AH7;AH8;AH9];
%还原对应工况下的压力/高度
P=U_P_ave+FI_P*AP;
H_l=U_H_ave+FI_H*AH;
%绘制图像
figure
plot(0:0.1:1199.9,P);
figure
plot(0:0.1:1199.9,H_l);