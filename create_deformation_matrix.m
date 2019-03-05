% segment roi region
I=imread('test.png');
[a, b, BW, xi, yi] = roipoly(I);

close all;

%%%
mask_pathname='.\Border_DB_Masks\';
[filename2, pathname2] = uigetfile({'*.jpg; *.tif; *.bmp; *.png'},'Load a mask image',mask_pathname);
mask=imread([pathname2 filename2]);

dim = size(BW);
[y x]=find(BW==1);
y=min(y);
imy=BW(y,:);
x=min(find(imy==1));

% detect boundaries and plot
boundary = bwtraceboundary(BW,[y, x],'N');%,'8','counterclockwise');  %x=column y=row
% boundary = bwtraceboundary(BW,[x, y],'N');
imshow(I);
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

 
 