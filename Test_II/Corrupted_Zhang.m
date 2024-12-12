function [output1,output2,output3] = Corrupted_Zhang(A0,n,k,z,p_value)
nk=n/k;B=zeros(k,k);
% Number of corruption nodes
zk=n*z;
% Corrupted adjacency matrix
A=A0;
A(nk-zk+1:nk,:)=A0(2*nk-zk+1:2*nk,:);
% Estemate the community connection probability matrix of the corrupted adjacency matrix
for i=1:k
B(i,i)=sum(sum(A((i-1)*nk+1:i*nk,(i-1)*nk+1:i*nk)))/(nk*(nk-1));
end
for i=2:k
    for j=1:i-1
B(i,j)=sum(sum(A((i-1)*nk+1:i*nk,(j-1)*nk+1:j*nk)))/(nk*nk);
    end
end
diag_B=diag(diag(B));B=B+B';B=B-diag_B;dist=1e-6; B=B+dist*ones(k,k);
% Estemate the entry-wise deviations of the corrupted adjacency matrix
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

Ln=max(max(abs(rho)));
Test_S=Ln^2-2*log(2*n*k)+log(log(2*n*k));
ev=-2*log(-2*sqrt(pi)*log(1-p_value));
if Test_S>ev
    Size_S=1;
    else
    Size_S=0;
end
output1=Size_S;
output2=Test_S;
output3=rho;

