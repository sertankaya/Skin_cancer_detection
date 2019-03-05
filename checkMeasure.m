function [ output_args ] = checkMeasure( fileName1,fileName2, title1,title2  )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%load matrix
close all;
%clc ;


Squ=importdata(fileName1);
Tri=importdata(fileName2);


[n,m]=size(Squ);
[n1,m1]=size(Tri);

I=eye(n,n);
I1=eye(n1,n1);

%%%%%%%%%%%%%%NORMALIZE THE MATRIX
Squ=Squ/norm(Squ);
Tri=Tri/norm(Tri);

disp('Matrix Norm 2:')
disp(title1);
disp(norm(Squ))
disp(title2);
disp(norm(Tri))


disp('Frobenious Norm:')
disp(title1);
N=Squ-I;
F=sqrt(trace(N*N'));

disp(F)


disp(title2);
N=Tri-I1;
F=sqrt(trace(N*N'));

disp(F);

disp(strcat(title1,'-Eigens'));
%[Squ_V,Squ_D]=eig(Squ);
%[Tri_V,Tri_D]=eig(Tri);
eigs(Squ)

disp(strcat(title2,'-Eigens'));
eigs(Tri)


 [U_Squ,S_Squ,V_Squ] = svd(Squ);
 [U_Tri,S_Tri,V_Tri] = svd(Tri);

 
 figure('name',strcat(title1,'-SVD'));
 plot(diag(S_Squ));
 figure('name',strcat(title1,'-Histogram'));
 hist(diag(S_Squ),n);
 
 figure('name',strcat(title2,'-SVD'));
 plot(diag(S_Tri));
 figure('name',strcat(title2,'-Histogram'));
 hist(diag(S_Tri),n);
 
 diag_S_Squ=diag(S_Squ);
 diag_S_Tri=diag(S_Tri);
 for i=2:n-2
   ddx_Squ(i)=(diag_S_Squ(i+1)-diag_S_Squ(i-1))/2.0;
   ddx_Tri(i)=(diag_S_Tri(i+1)-diag_S_Tri(i-1))/2.0;
   
 
 
 end
 
 
 figure('name',strcat(title1,'-_ddx'));
 plot(ddx_Squ);
 
 figure('name',strcat(title2,'-_ddx'));
 plot(ddx_Tri);
  
 figure('name',strcat('Laplacian_',title1));
 plot(abs(del2(diag(S_Squ-I))));
 disp(strcat('Laplacian std deivation for :',title1));
 %plot(abs(del2(diag(I))))
 std(del2(diag(S_Tri-I1)))
   
 figure('name',strcat('Laplacian_',title2));
 disp(strcat('Laplacian std deivation for :',title2));
 plot(abs(del2(diag(S_Tri-I1))))
 std(del2(diag(S_Squ-I)))
 
 
 startIndex=1;
 endIndex=n/2;
 totalIter=1;
%  
%  for k=1:totalIter
%      targetIndex=1;
%      sum=0;
%      for i=startIndex:endIndex
%        K2by2_Source=Squ(i:i+1,i:i+1)
%        K2by2_Target=Squ(endIndex+targetIndex:endIndex+targetIndex+1,endIndex+targetIndex:endIndex+targetIndex+1)
%        sum=sum+norm(K2by2_Source-K2by2_Target);
%        targetIndex=targetIndex+1;
%        
%     end
%  end
 

 

end

