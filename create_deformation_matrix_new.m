%------segment roi region
%%%%%  I=imread('test.png');
clc;clear all;close all

%images_pathname='.\Skin_Caner-100-Images\100_Images_JPG\';
images_pathname='C:\Users\Yasar\Desktop\Skin_Cancer-NN-100-Images\100_Images_JPG\';
[filename, pathname] = uigetfile({'*.jpg; *.tif; *.bmp; *.png'},'Load an ORIGINAL image',images_pathname);
img=imread([pathname filename]);
I=img;
%%%A=im2double(rgb2gray(img)); Because of red I commented out
% Red channel:
%A=im2double(img(:,:,3));

%%%% Mask
mask_pathname='C:\Users\Yasar\Desktop\Skin_Cancer-NN-100-Images\Border_DB_Masks\BD_DB_Masks';
[filename2, pathname2] = uigetfile({'*.jpg; *.tif; *.bmp; *.png'},'Load a MASK image',mask_pathname);
mask=imread([pathname2 filename2]);
mask=imresize(mask,[size(I,1) size(I,2)]);  % !!!! R C parametleri gir !!!!!!
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
%%%%%%%%%%%%%


dim = size(BW);
[y x]=find(BW==1);
y=min(y); 
imy=BW(y,:);
x=min(find(imy==1));

% detect boundaries and plot
boundary = bwtraceboundary(BW,[y, x],'N');%,'8','counterclockwise');  %x=column y=row
% boundary = bwtraceboundary(BW,[x, y],'N');
imshow(I);   %%% 4 satiri show icin geri ac
hold on;
plot(boundary(:,2),boundary(:,1),'g');%,'LineWidth',3);
hold on;

%contour = bwtraceboundary(BW, [row, col], 'W', 8, 50,...
%                                  'counterclockwise');

% find r of circle candidate
[k,t]=size(boundary);
r=k/(2*pi);

% Find the centroids of segmenteg region
s  = regionprops(BW, 'centroid');
centroids = cat(1, s.Centroid);
hold on;
plot(centroids(:,1), centroids(:,2), 'b+');
hold on;

%uzak=sqrt((centroids(:,2)-y)^2+(centroids(:,1)-x)^2);
%rc=centroids(:,2)-(r/uzak)*(centroids(:,2)-y)


% Draw circle using r and centroid points
th =linspace(0,2*pi,k);%0:pi/k:2*pi;
xunit = r * cos(th) + centroids(:,1);  % x axis corresponds to column      % [row col] 
yunit = r * sin(th) + centroids(:,2);  % y axis corresponds to row         % [row col]
h = plot(xunit, yunit,'r');
hold off

%unit=[yunit' xunit'];
% Create Deformation matrix
DG=zeros(k);
DG2=zeros(k);

for i=1:k
    DG(i,i)=sqrt( (boundary(i,1)-yunit(i))^2+(boundary(i,2)-xunit(i))^2);% euclidean
    DG2(i,i)=abs(boundary(i,1)-yunit(i))+abs(boundary(i,2)-xunit(i)); %manhattan
          % boundary(:,2)=colunm boundary(:,1)=row
end                                     %
%DG= pdist2(X,Y);
%DG= pdist2(boundary,unit);

H=[];
for i=1:k
    
    H=[H;boundary(i,1) boundary(i,2) DG(i,i) DG2(i,i)]; % [y=row x=column euclidean manhattan]
    
end

Excel = actxserver ('Excel.Application'); 
File='C:\Users\Yasar\Desktop\Gradient_Sinan_Abi_Kod\distance.xlsx'; 
if ~exist(File,'file') 
    ExcelWorkbook = Excel.workbooks.Add; 
    ExcelWorkbook.SaveAs(File,1); 
    ExcelWorkbook.Close(false); 
end 
%invoke(Excel.Workbooks,'Open',File);
ExcelWorkbook = Excel.workbooks.Open(File);

%filename = 'distance.xlsx';
format short
data = H;
xlswrite1(File,data,'Distance','B2');
ExcelWorkbook.Save 
ExcelWorkbook.Close(false) % Close Excel workbook. 
Excel.Quit; 
delete(Excel);





















%----------------------------------------------
% Place a circle at the centroid with the largest radius.
% boxX = centroids(:,1)- r;
% boxY = centroids(:,2)- r;
% 
% rectangle('Position',[ boxX boxY 2*r 2*r ],...
%     'Curvature',[1,1], 'EdgeColor', 'b', 'LineWidth', 3);
% hold on;

%col = round(dim(2)/2)-90;
%row = min(find(BW(:,col)));

 
 