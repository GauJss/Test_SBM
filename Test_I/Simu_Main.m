clear
% Model setting
n1=300; k0=2; k1=2;
n=n1*k0;
p0=0.5; p1=0.1; 
In_p=repelem(p0,k0);In_p=diag(In_p);
Out_p=repelem(p1,k0^2);Out_p=reshape(Out_p,k0,k0);Out_p=Out_p-diag(diag(Out_p));
B0=In_p+Out_p;
p_value=0.05;

loop=500; Boot_time=50; sample_time=100;    
% Empty matrix to save resluts 
time1=0;
Size_K=zeros(loop,1);Size_S=zeros(loop,1);Test_K=zeros(loop,1);Test_S=zeros(loop,1);
Size_KBoot=zeros(loop,1);Size_SBoot=zeros(loop,1);Test_KBoot=zeros(loop,1);Test_SBoot=zeros(loop,1);

Sizeaug_K=zeros(loop,1);Sizeaug_S=zeros(loop,1);Testaug_K=zeros(loop,1);Testaug_S=zeros(loop,1);
Sizeaug_KBoot=zeros(loop,1);Sizeaug_SBoot=zeros(loop,1);Testaug_KBoot=zeros(loop,1);Testaug_SBoot=zeros(loop,1);

P_Size_K=zeros(loop,1);P_Size_S=zeros(loop,1);P_Test_K=zeros(loop,1);P_Test_S=zeros(loop,1);
P_Size_KBoot=zeros(loop,1);P_Size_SBoot=zeros(loop,1); P_Test_KBoot=zeros(loop,1);P_Test_SBoot=zeros(loop,1);

P_Sizeaug_K=zeros(loop,1);P_Sizeaug_S=zeros(loop,1);P_Testaug_K=zeros(loop,1);P_Testaug_S=zeros(loop,1);
P_Sizeaug_KBoot=zeros(loop,1);P_Sizeaug_SBoot=zeros(loop,1); P_Testaug_KBoot=zeros(loop,1);P_Testaug_SBoot=zeros(loop,1);

position=zeros(loop,n);

parfor i=1:loop
tic;
disp(i) 
% Gerating adjacency matrix 
A=Generating(n,k0,B0);
% Network density
deg=sum(sum(A))/n^2;
sample_num=round(sqrt(1/deg)*(n/log(n))^(1/3));
% Use spectral k-means to estimate membership vector 
[hat_A,class]=SP_kmeans(A,k1);
position(i,:)=class;
% Estimate the community connection probability matrix and entry-wise deviations in Hu et al. (2021)
estype1='known';estype2='unknown';
[hat_SB,hat_Srho] = MLE(A,B0,n,k1,estype1,position(i,:));  % For test g
[hat_KB,hat_Krho] = MLE(A,B0,n,k1,estype2,position(i,:));  % For test K
% Testing results of the aug-bootstrap entry-wise deviations in Hu et al. (2021)
[Sizeaug_K(i,:),Sizeaug_S(i,:),Testaug_K(i,:),Testaug_S(i,:),hat_Krhoaug,hat_Srhoaug]=Zhang_Aug(A,hat_SB,n,k1,p_value);
[Sizeaug_KBoot(i,:),Sizeaug_SBoot(i,:),Testaug_KBoot(i,:),Testaug_SBoot(i,:)]=Zhang_AugBoot(Testaug_K(i,:),Testaug_S(i,:),hat_SB,n,k1,Boot_time,p_value); 
% Testing results of the sampling entry-wise deviations 
[P_Size_K(i,:),P_Size_S(i,:),P_Test_K(i,:),P_Test_S(i,:)] = Proposed(hat_Krho,hat_Srho,sample_time,sample_num,p_value);
[P_Size_KBoot(i,:),P_Size_SBoot(i,:),P_Test_KBoot(i,:),P_Test_SBoot(i,:)]=Proposed_Boot(P_Test_K(i,:),P_Test_S(i,:),hat_SB,n,k1,Boot_time,sample_time,sample_num,p_value);
% Testing results of the aug-bootstrap sampling entry-wise deviations 
[P_Sizeaug_K(i,:),P_Sizeaug_S(i,:),P_Testaug_K(i,:),P_Testaug_S(i,:)] = Proposed(hat_Krhoaug,hat_Srhoaug,sample_time,sample_num,p_value);
[P_Sizeaug_KBoot(i,:),P_Sizeaug_SBoot(i,:),P_Testaug_KBoot(i,:),P_Testaug_SBoot(i,:)]=Proposed_AugBoot(P_Testaug_K(i,:),P_Testaug_S(i,:),hat_SB,n,k1,Boot_time,sample_time,sample_num,p_value);

time1=toc;
end

AP_Pro=[sum(P_Size_K),sum(P_Size_KBoot),sum(P_Sizeaug_K),sum(P_Sizeaug_KBoot)]/loop; 



