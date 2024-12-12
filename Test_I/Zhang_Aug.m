function [output1,output2,output3,output4,output5,output6] = Zhang_Aug(A,B,n,k,p_value)
nk=n/k; n_new=n/k/2;
sizeB=length(B);
% calculate the edge probability of the newly added community
In_p=max(diag(B));
diag_B=diag(diag(B)); diag_B=reshape(diag_B,sizeB^2,1);
vec_B=reshape(B,sizeB^2,1);
Out_p=sort((vec_B-diag_B)/2); %% most of between probability are 0
Out_p=Out_p(find(Out_p>min(Out_p),1));
% Generate the newly added community
OA_new=zeros(n,n_new);
for i=1:k
    new_w1=binornd(1,Out_p,nk*n_new,1);
    OA_new((i-1)*nk+1:i*nk,1:n_new)=reshape(new_w1,nk,n_new);
end
new_w0=binornd(1,In_p,n_new*(n_new+1)/2,1); 
a1=length(new_w0);
a2=(-1+sqrt(1+8*a1))/2 ; 
new_w0a=zeros(a2,a2); 
for k1=1:a2
    for k2=k1:a2
        index=sum(a2:-1:a2-k1+2)+k2-k1+1;  
        new_w0a(k1,k2)=new_w0(index);      
    end
end
new_w0a=new_w0a+new_w0a'-diag(diag(new_w0a));
IA_new=new_w0a;
% Newyly adjacency Matrix
A_new=[A,OA_new;OA_new',IA_new];
A_new=A_new+diag(-diag(A_new));
% Estemate the community connection probability matrix of the new adjacency matrix
B_new=zeros(k+1,k+1);
for i=1:k
B_new(i,i)=sum(sum(A_new((i-1)*nk+1:i*nk,(i-1)*nk+1:i*nk)))/(nk*(nk-1));
end
B_new(k+1,k+1)=sum(sum(A_new(k*nk+1:k*nk+n_new,k*nk+1:k*nk+n_new)))/(n_new*(n_new-1));
for i=2:k
    for j=1:i-1
B_new(i,j)=sum(sum(A_new((i-1)*nk+1:i*nk,(j-1)*nk+1:j*nk)))/(nk*nk);
    end
end
for j=1:k
B_new(k+1,j)=sum(sum(A_new(k*nk+1:k*nk+n_new,(j-1)*nk+1:j*nk)))/(nk*n_new);
end
diag_B_new=diag(diag(B_new));
B_new=B_new+B_new';
B_new=B_new-diag_B_new;
dist=1e-6; B_new=B_new+dist*ones(k+1,k+1);
% Estemate the entry-wise deviations of the new adjacency matrix
S_rho_new=zeros(n+n_new,k+1);
for i=1:n+n_new
    i1=ceil(i/nk);
    for v=1:k
        if i1==v   
        S_rho_new(i,v)=sum((A_new(i,setdiff((v-1)*nk+1:v*nk,i))-B_new(i1,v))./sqrt(B_new(i1,v)*(1-B_new(i1,v))))/sqrt(nk-1);
        elseif i1<=k
        S_rho_new(i,v)=sum((A_new(i,(v-1)*nk+1:v*nk)-B_new(i1,v))./sqrt(B_new(i1,v)*(1-B_new(i1,v))))/sqrt(nk);
        else
        S_rho_new(i,k+1)=sum((A_new(i,setdiff(k*nk+1:k*nk+n_new,i))-B_new(i1,k+1))./sqrt(B_new(i1,k+1)*(1-B_new(i1,k+1))))/sqrt(n_new-1);
        end
    end  
end
% Use spectral k-means to estimate membership vector
[~,class]=SP_kmeans(A_new,k+1);  
position=class;
[~,K_rho_new] = MLE(A_new,B_new,n+n_new,k+1,'unknown',position); 
% Aug test Statistics in Hu et al. (2021)
S_Ln_new=max(max(abs(S_rho_new)));
K_Ln_new=max(max(abs(K_rho_new)));
S_Test_aug=S_Ln_new^2-2*log(2*(n+n_new)*(k+1))+log(log(2*(n+n_new)*(k+1)));
K_Test_aug=K_Ln_new^2-2*log(2*(n+n_new)*(k+1))+log(log(2*(n+n_new)*(k+1)));
ev=-2*log(-2*sqrt(pi)*log(1-p_value));
if S_Test_aug>ev
    S_Size_aug=1;
else
    S_Size_aug=0;
end
if K_Test_aug>ev
    K_Size_aug=1;
else
    K_Size_aug=0;
end
output1=K_Size_aug;
output2=S_Size_aug;
output3=K_Test_aug;
output4=S_Test_aug;
output5=K_rho_new;
output6=S_rho_new;
end

