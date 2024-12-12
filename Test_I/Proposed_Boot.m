function [output1,output2,output3,output4] = Proposed_Boot(TK,TS,B,n,k,m,sample_time,sample_num,p_value)
nk=n/k;
Test=zeros(m,1);
for loop=1:m
% Bootstrap procedure
A=Generating(n,k,B);
rho=zeros(n,k);
for i=1:n
    for v=1:k
        i1=ceil(i/nk);
        if i1==v
        rho(i,v)=sum((A(i,setdiff((v-1)*nk+1:v*nk,i))-B(i1,v))./sqrt(B(i1,v)*(1-B(i1,v))))/sqrt(nk-1);
        else
        rho(i,v)=sum((A(i,(v-1)*nk+1:v*nk)-B(i1,v))./sqrt(B(i1,v)*(1-B(i1,v))))/sqrt(nk);
        end
    end
end
% Calculate the sampling entry-wise deviations of the new adjacency matrix
samrho=zeros(sample_time,k);
for j=1:k
for i=1:sample_time
    sample=randsample(n,sample_num);
    sample_rho=rho(sample,j);
    samrho(i,j)=sum(sample_rho)/sqrt(sample_num);
end
end
%%
Ln=max(max(abs(samrho)));
Test(loop,:)=Ln^2-2*log(sample_time*k)+log(log(sample_time*k));
end
% MLE of the location parameters in the Gumble-distribution
param0=[-2*log(sqrt(pi)),2];
options = optimset('GradObj','off'); 
fun= {@(param)(m*log(param(2))-sum(-(Test-param(1))./param(2))+sum(exp(-(Test-param(1))./param(2))))};
[param_mle,~,~,~]=fminunc(fun,param0,options);
hat_mu=param_mle(1);hat_beta=param_mle(2);
%%
TK_Boot=-2*log(sqrt(pi))+2*(TK-hat_mu)/hat_beta; % For test K
TS_Boot=-2*log(sqrt(pi))+2*(TS-hat_mu)/hat_beta; % For test g
ev=-2*log(-sqrt(pi)*log(1-p_value));
if TK_Boot>ev
    Size_KB=1;
    else
    Size_KB=0;
end
if TS_Boot>ev
    Size_SB=1;
    else
    Size_SB=0;
end
output1=Size_KB;
output2=Size_SB;
output3=TK_Boot;
output4=TS_Boot;
end

