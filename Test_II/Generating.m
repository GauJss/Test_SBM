function output = Generating(n,k,B)  
A=zeros(n,n); nk=n/k; 
for i=1:k
    % Generate elements of the upper triangle of the adjacency matrix for
    % each block   
    w0=binornd(1,B(i,i),nk*(nk+1)/2,1); 
    a1=length(w0);
    % The dimension of the upper triangular matrix
    a2=(-1+sqrt(1+8*a1))/2 ; 
    % Transforms the vector into upper triangular
    w0a=zeros(a2,a2); 
for k1=1:a2
    for k2=k1:a2
        index=sum(a2:-1:a2-k1+2)+k2-k1+1;  
        w0a(k1,k2)=w0(index);      
    end
end
     w0a=w0a+w0a';
     w0a=w0a-diag(diag(w0a));
     A((i-1)*nk+1:i*nk,(i-1)*nk+1:i*nk)=w0a;
end
for i=2:k
    for j=1:i-1
    % Generate elements of the upper triangle of the adjacency matrix for
    % each block
    w1=binornd(1,B(i,j),nk*(nk+1)/2,1);
    a1=length(w1);
    % The dimension of the upper triangular matrix
    a2=(-1+sqrt(1+8*a1))/2 ; 
    % Transforms the vector into upper triangular
    w1a=zeros(a2,a2); 
for k1=1:a2
    for k2=k1:a2
        index=sum(a2:-1:a2-k1+2)+k2-k1+1;  
        w1a(k1,k2)=w1(index);      
    end
end
    w1a=w1a+w1a';
    w1a=w1a-diag(diag(w1a));
    A((i-1)*nk+1:i*nk,(j-1)*nk+1:j*nk)=w1a;  
    end
end   
% Symmetry
for i=1:n
    for j=1:i-1
        A(j,i)=A(i,j);
    end
end

A=A+diag(-diag(A));
output=A;

