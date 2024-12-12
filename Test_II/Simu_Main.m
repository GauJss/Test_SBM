clear
% Model setting
n1=300; k0=2; k1=2;
n=n1*k0;
z=0.05;
p0=0.5; p1=0.1; 
In_p=repelem(p0,k0);In_p=diag(In_p);
Out_p=repelem(p1,k0^2);Out_p=reshape(Out_p,k0,k0);Out_p=Out_p-diag(diag(Out_p));
B0=In_p+Out_p;
p_value=0.05;

loop=500; sample_time=100; Boot=50;
% Empty matrix to save resluts
time1=0;
Power_Size_S=zeros(loop,1);Power_P_Size_S=zeros(loop,1);Power_Test_S=zeros(loop,1);Power_P_Test_S=zeros(loop,1);
Power_Size_SBootAug=zeros(loop,1);Power_P_Size_SBootAug=zeros(loop,1);Power_Test_SBootAug=zeros(loop,1);Power_P_Test_SBootAug=zeros(loop,1);
position=zeros(loop,n);
parfor i=1:loop
tic;
disp(i) 
% Gerating adjacency matrix 
A=Generating(n,k0,B0);
% Network density
deg=sum(sum(A))/n^2;
sample_num=round(sqrt(1/deg)*(n/log(n))^(1/3));
% Testing results of the aug-bootstrap entry-wise deviations in Hu et al. (2021)
[Power_Size_S(i,:),Power_Test_S(i,:),rho]=Corrupted_Zhang(A,n,k1,z,p_value);
[Power_Size_SBootAug(i,:),Power_Test_SBootAug(i,:),rho_BootAug,hatB_BootAug]=Corrupted_Zhang_BootAug(A,n,k1,z,Boot,p_value);
% Testing results of the aug-bootstrap sampling entry-wise deviations 
[Power_P_Size_S(i,:),Power_P_Test_S(i,:)]=Corrupted_Proposed(rho,sample_time,sample_num,p_value);
[Power_P_Size_SBootAug(i,:),Power_P_Test_SBootAug(i,:)]= Corrupted_Proposed_BootAug(rho_BootAug,hatB_BootAug,sample_time,sample_num,n,Boot,p_value);
time1=toc;
end

a=[sum(Power_P_Size_S),sum(Power_P_Size_SBootAug),sum(Power_Size_S),sum(Power_Size_SBootAug)]/loop;

