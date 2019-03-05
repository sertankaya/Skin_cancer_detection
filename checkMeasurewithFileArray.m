clear all;
close all;
clc;

fileList={
'6094d.txtdeformation.csv',
'Ma0899ad.txtdeformation.csv',
'Mm0899ad.txtdeformation.csv',
'jb0298ad.txtdeformation.csv',
'mj0798ad.txtdeformation.csv',
'4056d.txtdeformation.csv',
'4058d.txtdeformation.csv',
'As0598ad.txtdeformation.csv',
'zr0798ad.txtdeformation.csv',
'6714d.txtdeformation.csv',
'6688d.txtdeformation.csv',
'4063d.txtdeformation.csv',
'6599d.txtdeformation.csv',
'6089d.txtdeformation.csv',
'6031d.txtdeformation.csv',
'6005d.txtdeformation.csv',
'6091d.txtdeformation.csv',
'Wa0499ad.txtdeformation.csv',
'6077d.txtdeformation.csv',
'4023.txtdeformation.csv',
'6642d.txtdeformation.csv',
'kh0598.txtdeformation.csv',
'6050d.txtdeformation.csv',
'6677d.txtdeformation.csv',
'hm1096ad.txtdeformation.csv',
'6687d.txtdeformation.csv',
'6075d.txtdeformation.csv',
'Ld1196ad.txtdeformation.csv',
'6644d.txtdeformation.csv',
'6074d.txtdeformation.csv',
'6028d.txtdeformation.csv',
'6605d.txtdeformation.csv',
'6084d.txtdeformation.csv',
'6023d.txtdeformation.csv',
'6602d.txtdeformation.csv',
'6044d.txtdeformation.csv',
'6709d.txtdeformation.csv',
'6612d.txtdeformation.csv',
'6051d.txtdeformation.csv',
'6025d.txtdeformation.csv',
'6024d.txtdeformation.csv',
'6017d.txtdeformation.csv',
'6079d.txtdeformation.csv',
'6643d.txtdeformation.csv',
'Kr1196ad.txtdeformation.csv',
'6622d.txtdeformation.csv',
'4025d.txtdeformation.csv',
'6706d.txtdeformation.csv',
'6054d.txtdeformation.csv',
'6102d.txtdeformation.csv',
'6055d.txtdeformation.csv',
'6049d.txtdeformation.csv',
'6637d.txtdeformation.csv',
'6093d.txtdeformation.csv',
'6045d.txtdeformation.csv',
'6707d.txtdeformation.csv',
'Fj0297.txtdeformation.csv',
'6053d.txtdeformation.csv',
'6627d.txtdeformation.csv',
'6061d.txtdeformation.csv',
'6065d.txtdeformation.csv',
'6617d.txtdeformation.csv',
'6675d.txtdeformation.csv',
'6667d.txtdeformation.csv'};

[n,m]=size(fileList);
%checkMeasure('Star.csv','Arrow.csv','STAR','ARROW');
for i=1:n/2 
  x=(i-1)*2+1;
  y=(i-1)*2+2;
  if(y>n)y=n;end
  fileName1=fileList{x};
  fileName2=fileList{y};
 

  fprintf('\n');
  checkMeasure(fileName1,fileName2,fileName1,fileName2);
  input('Press a key to continue');
end

 


