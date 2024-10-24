clc

[FI_P,A_P,LAM_P,E_P,n_P]=snapshotPOD(U_P_wavy);
energy(E_P,n_P);


[FI_H,A_H,LAM_H,E_H,n_H]=snapshotPOD(U_H_wavy);
energy(E_H,n_H);

%降阶后，手动复制结果