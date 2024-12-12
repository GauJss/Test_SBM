function [output1,output2] = SP_kmeans(A,k)
% Eigenvalue decomposition
[eig_vector,eigvalue]=eig(A);
% The k largest eigenvalues and their corresponding eigenvectors
[~,order]=sort(diag(eigvalue),'descend'); 
eig_Lead = eig_vector(:,order);
eig_Lead=eig_Lead(:,1:k); 
% The estimated membership vector
hat_sigma = kmeans(eig_Lead,k);
num=0;
for kk=1:k
    num1=find(hat_sigma==kk);
    num=[num;num1];
end
% The newly sorted adjacency matrix
num(1)=[];
hat_A=A(num,num);
output1=hat_A;
output2=hat_sigma;
end

