% Original
clc;clear all;
close all;

images_pathname='.\';
[filename, pathname] = uigetfile({'*.jpg; *.tif; *.bmp; *.png'},'Load an ORIGINAL image',images_pathname);
img=imread([pathname filename]);
A=im2double(rgb2gray(img));


%%%% Mask
mask_pathname='.\';
[filename2, pathname2] = uigetfile({'*.jpg; *.tif; *.bmp; *.png'},'Load a MASK image',mask_pathname);
mask=imread([pathname2 filename2]);
mask=imresize(mask,[size(A,1) size(A,2)]);  % !!!! R C parametleri gir !!!!!!
BW=imcomplement(mask);
BW=im2bw(BW);
%BW=bwareaopen(BW, 8700);
BW=bwareaopen(BW, 13000);

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
i=1;
j=1;
flag=1;
homogenty=[];
r=10;
ang=0:0.01:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);
figure;imshow(img);
hold on;
plot(boundary(:,2),boundary(:,1),'g');%,'LineWidth',3);
hold on;
%close all;
while (flag==1)

plot(boundary(j,2)+xp,boundary(i,1)+yp);
%hold off;
%close all
%%%  end

tmp=zeros(size(A,1),size(A,2));
figure;imshow(tmp);hold on;
plot(boundary(j,2)+xp,boundary(i,1)+yp,'w');
%roifill(I,boundary(j,2)+xp(1),boundary(i,1)+yp(1));
% roifill(new,boundary(j,2),boundary(i,1));
%J = roifill(I,c,r);
u=getframe;
new=im2bw(u.cdata);
%new=imresize(new,[size(A,1),size(A,2)]);
se = strel('disk',10);
BW = imclose(new,se);
BW=imresize(BW,[size(A,1),size(A,2)]);
imshow(BW);
%close figure 2;



%  clear outside
percentage_of_edges_to_clip_top=5;
percentage_of_edges_to_clip=10;
[h,w]=size(BW);
left_mask_boundary=1+(w/2*percentage_of_edges_to_clip/100);
right_mask_boundary=w-(w/2*percentage_of_edges_to_clip/100);
top_mask_boundary=1+(h/2*percentage_of_edges_to_clip_top/100);
bottom_mask_boundary=h-(h/2*percentage_of_edges_to_clip/100);

BW(:,1:left_mask_boundary)=0;
BW(:,right_mask_boundary:end)=0;
BW(1:top_mask_boundary,:)=0;
BW(bottom_mask_boundary:end,:)=0;
imshow(BW);


result=A.*BW;
imshow(result);
vector=A(find(BW))';

%stats = graycoprops(result,'homogeneity'); % Homoggoneity for C
ans = graycoprops(graycomatrix(vector), 'Homogeneity');
homogenty=[homogenty ans.Homogeneity];

i=round(i+(3*r/2));
j=round(j+(3*r/2));

   if (i>=size(boundary,1)|| j>=size(boundary,1))
     flag=2;
   end
   close figure 2
   %close figure 3
   hold on;
end



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