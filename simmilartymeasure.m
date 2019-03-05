%load matrix
clear all;
close all;
clc;

%Squ=importdata('Square.csv');
%Tri=importdata('Triangle.csv');


Squ=importdata('Star.csv');
Tri=importdata('Arrow.csv');

[n,m]=size(Tri);

I=eye(n,n);

disp('Matrix Norm 2')
disp('Square ');
disp(norm(Squ))
disp('Triangle');
disp(norm(Tri))


disp('Frobenious Norm')
disp('Square ');
N=Squ-I;
F=sqrt(trace(N*N'));

disp(F)


disp('Triangle ');
N=Tri-I;
F=sqrt(trace(N*N'));

disp(F);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tansel Halic
% I defined the similiary measure as an eigen value problem. The derivation
% of the similiarity problem is here below
% We might need to look into the eigen values and vectors of deformation
% gradient marix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Derivation%%%%%%%%%%%%%%%%%%%%%%
%Two matrices A and B are called to be similar if there exists a transformation defined as the following;  
%A=inver of S *B*S
%In our case, we define the similarity measure with respect to the closeness  to the perfect similarity. 
%Therefore, we’re looking a deformation gradient that is close to the perfect symmetry that is a  ? dimension circle.
%Therefore our problem  becomes;
% ?I=inverse of S *B*S
% ?S= BS
% BS- ?S=0
% (B- ?I)S=0
% B- ?I=0
%This will become (B- ?I)v=0. Our problem simply becomes Eigen value problem.
%For non-trivial solution , the det(B- ?I) needs to be zero. 

%This statement below needs to be validated as well.
%Therefore, the similarity measure also confirms that eigen values do exist and the
%distributions of the eigen values with respect to the ? gives use a similarity measure.  


%%basic check of eigen values of difference deformation gradient.
%[Squ_V,Squ_D]=eig(Squ);
%[Tri_V,Tri_D]=eig(Tri);
%eig(Squ)


%%idea of the SVD check. more info is in the paper. we need to look into
%%this more by comparing the angles as well
 [U_Squ,S_Squ,V_Squ] = svd(Squ);
 [U_Tri,S_Tri,V_Tri] = svd(Tri);
 figure('name','Squ-SVD');
 plot(diag(S_Squ));
 figure('name','Squ-Hist');
 hist(diag(S_Squ),n);
 
 figure('name','Tri-SVD');
 plot(diag(S_Tri));
 figure('name','Tri-Hist');
 hist(diag(S_Tri),n);
 
 diag_S_Squ=diag(S_Squ);
 diag_S_Tri=diag(S_Tri);
 for i=2:n-2
   ddx_Squ(i)=(diag_S_Squ(i+1)-diag_S_Squ(i-1))/2.0;
   ddx_Tri(i)=(diag_S_Tri(i+1)-diag_S_Tri(i-1))/2.0;
   
 
 
 end
 
 
 figure('name','Squ_ddx');
 plot(ddx_Squ);
 
 figure('name','Tri-ddx');
 plot(ddx_Tri);
  
 figure('name','Squ_Laplacian');
 plot(abs(del2(diag(S_Squ-I))));
 disp('Laplacian std deivation for Square ');
 %plot(abs(del2(diag(I))))
 std(del2(diag(S_Tri-I)))
   
 figure('name','Tri-Laplacian');
 disp('Laplacian std deivation for Tri ');
 plot(abs(del2(diag(S_Tri-I))))
 std(del2(diag(S_Squ-I)))
 
 
 
 
 
 
 
 