U_P_wavy=U_P-repmat(U_P_ave,1,9900);
ERROR=[];
for L=200:200
    for N=1:9900
        error=norm(U_P_wavy-FI_P(:,1:L)*(FI_P(:,1:L).')*U_P_wavy)/norm(U_P(:,N));
        ERROR=[ERROR;error];
        [L,N]
    end
end
error_ave=mean(ERROR);
