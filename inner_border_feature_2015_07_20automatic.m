% Original
clc;clear all;
close all;
format short

d1=size(dir('C:\Users\Sertan\Desktop\Skin_Cancer-NN-100-Images\100_Images_JPG\'),1);

Excel = actxserver ('Excel.Application'); 
File='C:\Users\Sertan\Desktop\Gradient_Sinan_Abi_Kod\data.xlsx'; 
if ~exist(File,'file') 
    ExcelWorkbook = Excel.workbooks.Add; 
    ExcelWorkbook.SaveAs(File,1); 
    ExcelWorkbook.Close(false); 
end 
%invoke(Excel.Workbooks,'Open',File);
ExcelWorkbook = Excel.workbooks.Open(File);


for index=6:d1  %d1 e - 3 ten basla cunku 1, 2 (. ve ..)

ip='C:\Users\Sertan\Desktop\Skin_Cancer-NN-100-Images\100_Images_JPG\';
%images_pathname='.\Skin_Caner-100-Images\100_Images_JPG\';
images_pathname=dir('C:\Users\Sertan\Desktop\Skin_Cancer-NN-100-Images\100_Images_JPG\');
%[filename, pathname] = uigetfile({'*.jpg; *.tif; *.bmp; *.png'},'Load an ORIGINAL image',images_pathname);

img=imread([ip images_pathname(index).name]);
A=im2double(rgb2gray(img)); %Because of red I commented out
% Red channel:
%%% YCBCR
%YCBCR = rgb2ycbcr(img); %%% bitmedi bak buna !!
%A=im2double(YCBCR(:,:,3));

%%% HSV
%hsv=rgb2hsv(img); 
%A=im2double(hsv(:,:,3));



%%%% Mask
mp='C:\Users\Sertan\Desktop\Skin_Cancer-NN-100-Images\Border_DB_Masks\';
masks_pathname=dir('C:\Users\Sertan\Desktop\Skin_Cancer-NN-100-Images\Border_DB_Masks');

mask=imread([mp masks_pathname(index).name]);
mask=imresize(mask,[size(A,1) size(A,2)]);  % !!!! R C parametleri gir !!!!!!
BW=imcomplement(mask);
BW=im2bw(BW);

%%%%%%%% CLEAR Outside of the region
%  clear outside
percentage_of_edges_to_clip_top=10;
percentage_of_edges_to_clip=10;
[h,w]=size(BW);
left_mask_boundary=floor(1+(w/2*percentage_of_edges_to_clip/100));
right_mask_boundary=floor(w-(w/2*percentage_of_edges_to_clip/100));
top_mask_boundary=floor(1+(h/2*percentage_of_edges_to_clip_top/100));
bottom_mask_boundary=floor(h-(h/2*percentage_of_edges_to_clip_top/100)); %% Bu kisim onemli, eger sekil at ve ustte yakinsa bu degismeli
BW(:,1:left_mask_boundary)=0;
BW(:,right_mask_boundary:end)=0;
BW(1:top_mask_boundary,:)=0;
BW(bottom_mask_boundary:end,:)=0;
%%%%%%%%%%%%%%%

%BW=bwareaopen(BW, 8700);   %%%% Burda da 2 secenek olmali, image gore
% tmp=round(h*w/21);  % gui skindeki degeri
tmp=round(h*w/33);
BW=bwareaopen(BW, tmp);
%BW=bwareaopen(BW, 18000);



dim = size(BW);
[y x]=find(BW==1);
y=min(y);
imy=BW(y,:);
x=min(find(imy==1));

% detect boundaries and plot
boundary = bwtraceboundary(BW,[y, x],'N');%,'8','counterclockwise');  %x=column y=row
% boundary = bwtraceboundary(BW,[x, y],'N');
% imshow(img);
% hold on;
% plot(boundary(:,2),boundary(:,1),'g');%,'LineWidth',3);
% hold on;


%%%%%%%%%%%% DRAW CIRCLE :

%%%   function circle(x,y,r)
%x and y are the coordinates of the center of the circle
%r is the radius of the circle
%0.01 is the angle step, bigger values will draw the circle faster but
%you might notice imperfections (not very smooth)

flag=1;

y=[5 7 10 15];
%alphabet = {'A','G','M','S'};
alphabet = {'B','H','N','T'};

for s=1:4
homogenty=[];
ort=[];
sapma=[];
border_point=[];
nd=[];
u=[];
q=[];
w=[];
mask2=[];
i=1;
j=1;
%%%%%%%%%%%%%%%%%%%
r=y(s);      %%%% mavi
ang= linspace(0,2*pi,size(boundary,1));
xp=r*cos(ang);
yp=r*sin(ang);
%%%%%%%%%%%%%%%%
%%% figure;imshow(img); because of red I commented out
%%%%%%% Red channel
figure;imshow(A);
hold on;
plot(boundary(:,2),boundary(:,1),'g');%,'LineWidth',3);
hold on;
%%%%%%%%%%%%%
% Find the centroids of segmenteg region
x  = regionprops(BW, 'centroid');
centroids = cat(1, x.Centroid);
hold on;
plot(centroids(:,1), centroids(:,2), 'b+');
hold on;

%Create a logical image of a circle with specified
% diameter, center, and image size.
% First create the image.

%[columnsInImage rowsInImage] = meshgrid(1:size(A,2), 1:size(A,1));


while (flag==1)
    
%tmp=[boundary(j,1) boundary(j,2)];
border_point=[boundary(j,2) boundary(j,1)];  % boundary(j,2)=X boundary(i,1)=Y

% find itself divide by norm
nd=(([centroids-border_point])./norm([centroids-border_point])); % division element by element

% multiply (element by element) it with boundary(X,Y)coordinates and multiply(element by element) the result with r
u=[boundary(j,2)+nd(1,1).*r boundary(j,1)+nd(1,2).*r]; % new coordinates (X,Y)

%%%%%%%%%%%%%%%%%%%%%%
    
%%% CIRCLE CIZIMI GEREKIYOR
%plot(boundary(j,2)+xp,boundary(i,1)+yp); %%%Mavi
plot(u(1,1)+xp,u(1,2)+yp); %%%Mavi

%%%%%%%%%%%%%%%%%

% Next create the circle in the image.
%centerX = boundary(j,2);
%centerY = boundary(i,1);
centerX = u(1,1);
centerY = u(1,2);

%%%%% Polygon codes

x2=u(1,1)+xp;
y2=u(1,2)+yp;
[q,w] = polybool('intersection',boundary(:,2),boundary(:,1),x2',y2');

%[row, col] = find(isnan(k));
q = fill_nans(q);
w = fill_nans(w);


mask2 = poly2mask(q,w,size(A,1),size(A,2));
 
plot(q, w, 'r');

figure;imshow(mask2)

%polybool()
radius =r;
r=radius;
BWn=mask2;
%figure;imshow(BWn);

%  clear outside
percentage_of_edges_to_clip_top=5;
percentage_of_edges_to_clip=10;
[h,w]=size(BWn);
left_mask_boundary=floor(1+(w/2*percentage_of_edges_to_clip/100));
right_mask_boundary=floor(w-(w/2*percentage_of_edges_to_clip/100));
top_mask_boundary=floor(1+(h/2*percentage_of_edges_to_clip_top/100));
bottom_mask_boundary=floor(h-(h/2*percentage_of_edges_to_clip_top/100)); %% Bu kisim onemli, eger sekil at ve ustte yakinsa bu degismeli

BWn(:,1:left_mask_boundary)=0;
BWn(:,right_mask_boundary:end)=0;
BWn(1:top_mask_boundary,:)=0;
BWn(bottom_mask_boundary:end,:)=0;
imshow(BWn);


result=A.*BWn;
imshow(result);
vector=A(find(BWn))';

%stats = graycoprops(result,'homogeneity'); % Homoggoneity for C
ans = graycoprops(graycomatrix(vector), 'Homogeneity');
homogenty=[homogenty ans.Homogeneity];

ort=[ort mean(vector)];
sapma=[sapma std(vector)];


i=round(i+(3*r/2));
j=round(j+(3*r/2));

   if (i>=size(boundary,1)|| j>=size(boundary,1))
       %if s>4
          flag=2;
       %end
   end
   close figure 2
   %close figure 3
   hold on; %%Mavi
end

if s<4
    flag=1;
elseif s==4;
    flag=2;
 else
    flag=2;
 end

clc
%display(['Average Homogeneity ' num2str(ycap)])
homo = fill_nans(homogenty');
orto = fill_nans(ort');
sapo = fill_nans(sapma');

homogenty=homo';
ort=orto';
sapma=sapo';



display('Average Homogeneity')
mean(homogenty)
display('Minumum Homogeneity')
min(homogenty)

display(' ');
display('Average Mean')
mean(ort)
display('Minumum Mean ')
min(ort)

display(' ');
display('Average Std')
mean(sapma)
display('Minumum Std ')
min(sapma)


filename = 'data.xlsx';
format short
data = [mean(homogenty) min(homogenty) mean(ort) min(ort) mean(sapma) min(sapma)];
xlRange =strcat(alphabet(s), num2str(index-2));
xlRange=num2str(cell2mat(xlRange));
%location=[alphabet(s) num2str(index-2)];

xlswrite1(File,data,'Sheet1',xlRange);
%%%xlswrite1(File,data,location);
ExcelWorkbook.Save 


close figure 1
end
onm =strcat('A', num2str(index-2));
%onm=num2str(cell2mat(onm));
isim={images_pathname(index).name};
xlswrite1(File,isim,'Sheet1',onm);


end

ExcelWorkbook.Save 
ExcelWorkbook.Close(false) % Close Excel workbook. 
Excel.Quit; 
delete(Excel);
















% %%%%%%%%%%%%%%
% %%%%%%%%%%%%%
% % find r of circle candidate
% %[k,t]=size(boundary);
% %r=k/(2*pi);
% r=10;
% % Find the centroids of segmenteg region
% s  = regionprops(BW, 'centroid');
% centroids = cat(1, s.Centroid);
% hold on;
% plot(centroids(:,1), centroids(:,2), 'b+');
% hold on;
% 
% %uzak=sqrt((centroids(:,2)-y)^2+(centroids(:,1)-x)^2);
% %rc=centroids(:,2)-(r/uzak)*(centroids(:,2)-y)
% 
% 
% % Draw circle using r and centroid points
% th =linspace(0,2*pi,k);%0:pi/k:2*pi;
% xunit = r * cos(th) + centroids(:,1);  % x axis corresponds to column      % [row col] 
% yunit = r * sin(th) + centroids(:,2);  % y axis corresponds to row         % [row col]
% h = plot(xunit, yunit,'r');
% hold off

% f=getframe(gcf);
% new=f.cdata;
% imwrite(new,'new.jpg');

%roifill(I,boundary(j,2)+xp(1),boundary(i,1)+yp(1));
% roifill(new,boundary(j,2),boundary(i,1));
%J = roifill(I,c,r);