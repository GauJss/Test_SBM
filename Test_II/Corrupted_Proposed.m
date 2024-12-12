function [output1,output2] = Corrupted_Proposed(rho,sample_time,sample_num,p_value)
sizerho=size(rho);
n=sizerho(1);k=sizerho(2);
% Calculate the sampling entry-wise deviations
Srho=zeros(sample_time,k);
for j=1:k
for i=1:sample_time
    sample=randsample(n,sample_num);
    sample_rho=rho(sample,j);
    Srho(i,j)=sum(sample_rho)/sqrt(sample_num);
end
end
% Compare with critical value
Ln=max(max(abs(Srho)));
Test_S=Ln^2-2*log(sample_time*k)+log(log(sample_time*k));
ev=-2*log(-sqrt(pi)*log(1-p_value));
if Test_S>ev
    Size_S=1;
    else
    Size_S=0;
end
output1=Size_S;
output2=Test_S;
end
