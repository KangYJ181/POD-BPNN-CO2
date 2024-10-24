function energy(E,p)
% 功能:绘制POD基相对能量分布图
% E:POD基相对能量列向量
% p:前p个POD基
yyaxis left % 画Ei分布图
bar(E(1:p),'red');
xlim([0.5,p+0.5]);
ylim([0,1.1]);
xlabel('Number of POD basis');
ylabel('Relative energy');
yyaxis right % 画∑Ei分布图
E_sum=[];
for k=1:p
    E_sum=[E_sum,sum(E(1:k))];
end
plot(1:p,E_sum,'-bs','MarkerFaceColor','b','LineWidth',1.5);
ylim([0,1.1]);
ylabel('Accumulated relative energy');
end