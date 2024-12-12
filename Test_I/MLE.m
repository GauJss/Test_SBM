function [output1,output2] = MLE(A,B0,n,k,estype,position)
nk=n/k;B=zeros(k,k);
% Estimate community connection probability matrix 
if strcmpi(estype , 'known') % true membership vector g is known
for i=1:k
B(i,i)=sum(sum(A((i-1)*nk+1:i*nk,(i-1)*nk+1:i*nk)))/(nk*(nk-1));
end
for i=2:k
    for j=1:i-1
B(i,j)=sum(sum(A((i-1)*nk+1:i*nk,(j-1)*nk+1:j*nk)))/(nk*nk);
    end
end
elseif strcmpi(estype , 'unknown') % true membership vector g is unknown 
for i=1:k
pos=find(position==i);
B(i,i)=sum(sum(A(pos,pos)))/(length(pos)*(length(pos)-1));
end
for i=2:k
    for j=1:i-1
        pos1=find(position==i); pos2=find(position==j);
        B(i,j)=sum(sum(A(pos1,pos2)))/(length(pos1)*length(pos2));
    end
end
else
    error('The case does not exit!');    
end

% symmetry
diag_B=diag(diag(B));
B=B+B';
B=B-diag_B;
dist=1e-6; B=B+dist*ones(k,k);

%% Calculate entry-wise deviations
rho=zeros(n,k);
if strcmpi(estype , 'known') % true membership vector g is known 
  for i=1:n
      for v=1:k
        i1=ceil(i/nk);
        if i1==v
        rho(i,v)=sum((A(i,setdiff((v-1)*nk+1:v*nk,i))-B0(i1,v))./sqrt(B0(i1,v)*(1-B0(i1,v))))/sqrt(nk-1);
        else
        rho(i,v)=sum((A(i,(v-1)*nk+1:v*nk)-B0(i1,v))./sqrt(B0(i1,v)*(1-B0(i1,v))))/sqrt(nk);
        end
      end
  end
elseif strcmpi(estype , 'unknown') % true membership vector g is unknown 
    for i=1:n
    for v=1:k
        i1=position(i);
        pos=find(position==v);
        if i1==v 
         rho(i,v)=sum((A(i,setdiff(pos,i))-B(i1,v))./sqrt(B(i1,v)*(1-B(i1,v))))/sqrt(length(pos)-1);
        else
         rho(i,v)=sum((A(i,setdiff(pos,i))-B(i1,v))./sqrt(B(i1,v)*(1-B(i1,v))))/sqrt(length(pos));
        end
    end
    end
else
    error('The case does not exit!'); 
end
%%
output1=B;
output2=rho;


