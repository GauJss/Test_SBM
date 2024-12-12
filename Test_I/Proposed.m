function [output1,output2,output3,output4] = Proposed(hat_Krho,hat_Srho,sample_time,sample_num,p_value)
sizerho=size(hat_Krho);
n=sizerho(1);k=sizerho(2);
% Calculate the sampling entry-wise deviations of the new adjacency matrix
krho=zeros(sample_time,k);srho=zeros(sample_time,k);
for j=1:k
for i=1:sample_time
    sample=randsample(n,sample_num);
    sample_krho=hat_Krho(sample,j);
    sample_srho=hat_Srho(sample,j);
    krho(i,j)=sum(sample_krho)/sqrt(sample_num); % For test K
    srho(i,j)=sum(sample_srho)/sqrt(sample_num); % For test g
end
end

Ln_K=max(max(abs(krho)));
Test_K=Ln_K^2-2*log(sample_time*k)+log(log(sample_time*k));
ev=-2*log(-sqrt(pi)*log(1-p_value));
if Test_K>ev
    Size_K=1;
else
    Size_K=0;
end

Ln_S=max(max(abs(srho)));
Test_S=Ln_S^2-2*log(sample_time*k)+log(log(sample_time*k));
ev=-2*log(-sqrt(pi)*log(1-p_value));
if Test_S>ev
    Size_S=1;
    else
    Size_S=0;
end
output1=Size_K;
output2=Size_S;
output3=Test_K;
output4=Test_S;
end
