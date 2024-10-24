function [FI,A,LAM,E,n]=snapshotPOD(U_NP)
% 功能:实现快照POD,适用于变量列数P<<不变量行数N的矩阵(例如矩阵行指标对应空间坐标，列指标对应工况)
% U_NP:数据矩阵,行数N,列数P
% LAM:按绝对值降序排列的特征值列向量
% FI:根据LAM排列的特征向量矩阵FI=[FI1,FI2,...]
% A:POD基系数
% E:模态相对能量列向量
% n:模态个数
C=(U_NP.')*U_NP;
[FI,LAM]=eig(C);
FI=U_NP*FI;
[LAM,index]=sort(diag(abs(LAM)),'descend');                                % 将特征值按绝对值降序重排
FI=FI(:,index);                                                            % 重排C的特征向量矩阵
[~,n]=size(FI);                                                            % n:POD基的个数
for k=1:n
    FI(:,k)=FI(:,k)/norm(FI(:,k));                                         % 将所有POD基归一化
end
A=(FI.')*U_NP;                                                             % 求POD基系数矩阵
E=LAM/sum(LAM);                                                            % 求各模态的相对能量
end