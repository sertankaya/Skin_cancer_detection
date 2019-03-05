%------segment roi region
%%%%%  I=imread('test.png');
clc;clear all;close all

%images_pathname='.\Skin_Caner-100-Images\100_Images_JPG\';
images_pathname='C:\Users\Sertan\Desktop\Skin_Cancer-NN-100-Images\100_Images_JPG\';
[filename, pathname] = uigetfile({'*.jpg; *.tif; *.bmp; *.png'},'Load an ORIGINAL image',images_pathname);
img=imread([pathname filename]);
I=img;
%%%A=im2double(rgb2gray(img)); Because of red I commented out
% Red channel:
%A=im2double(img(:,:,3));

%%%% Mask
mask_pathname='C:\Users\Sertan\Desktop\Skin_Cancer-NN-100-Images\Border_DB_Masks\BD_DB_Masks';
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
% plot(boundary(i,2),boundary(i,1),'+r');
% plot(boundary(i+4,2),boundary(i+4,1),'*b');
angle41=[];
angle42=[];
angle43=[];
angle44=[];

for i=1:12:size(boundary,1)-12
    angle41=[angle41; atand(((boundary(i+12,1)-boundary(i,1))/(boundary(i+12,2)-boundary(i,2))))]; 
end

for j=2:12:size(boundary,1)-12
    angle42=[angle42; atand(((boundary(j+12,1)-boundary(j,1))/(boundary(j+12,2)-boundary(j,2))))]; 
end

for j=3:12:size(boundary,1)-12
    angle43=[angle43; atand(((boundary(j+12,1)-boundary(j,1))/(boundary(j+12,2)-boundary(j,2))))]; 
end

for j=4:12:size(boundary,1)-12
    angle44=[angle44; atand(((boundary(j+12,1)-boundary(j,1))/(boundary(j+12,2)-boundary(j,2))))]; 
end

close all;
t1=1:size(angle41',2);
t2=1:size(angle42',2);
t3=1:size(angle43',2);
t4=1:size(angle44',2);

display('Std ler');
display('Std Angle121');
std(angle41)
display('Std Angle122');
std(angle42)
display('Std Angle123');
std(angle43)
display('Std Angle124');
std(angle44)

display('Std lerin Mean');
(std(angle41)+std(angle42)+std(angle43)+std(angle44))/4


figure
plot(t1,angle41','--rs',...
    'LineWidth',1.5,...
    'MarkerSize',5,...
    'MarkerEdgeColor','g',...
    'MarkerFaceColor',[0.3,0.3,0.3])
title(['Angles  ' filename]);
xlabel('Distance 4 pixels apart');
ylabel('Angles');


figure
plot(t2,angle42','--ks',...
    'LineWidth',1.5,...
    'MarkerSize',5,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5])
title(['Absolute Value Angles Averages ' filename]);
xlabel('Distance 4 pixels apart');
ylabel('Angles');

figure
plot(t3,angle43','--bs',...
    'LineWidth',1.5,...
    'MarkerSize',5,...
    'MarkerEdgeColor','m',...
    'MarkerFaceColor',[0.5,0.5,0.5])
title(['Angles' filename]);
xlabel('Distance 4 pixels apart');
ylabel('Angles ');

figure
plot(t4,angle44','--bs',...
    'LineWidth',1.5,...
    'MarkerSize',5,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor',[0.5,0.5,0.5])
title(['Angles' filename]);
xlabel('Distance 4 pixels apart');
ylabel('Angles ');      

















%%%%%%%%%%%%%%%%%% 4 pixels apart %%%%%%%%%%%%%%%%%%%%%5


% % for i=1:4:size(boundary,1)-4
% %     angle41=[angle41; atand(((boundary(i+4,1)-boundary(i,1))/(boundary(i+4,2)-boundary(i,2))))]; 
% % end
% % 
% % for j=2:4:size(boundary,1)-4
% %     angle42=[angle42; atand(((boundary(j+4,1)-boundary(j,1))/(boundary(j+4,2)-boundary(j,2))))]; 
% % end
% % 
% % for j=3:4:size(boundary,1)-4
% %     angle43=[angle43; atand(((boundary(j+4,1)-boundary(j,1))/(boundary(j+4,2)-boundary(j,2))))]; 
% % end
% % 
% % for j=4:4:size(boundary,1)-4
% %     angle44=[angle44; atand(((boundary(j+4,1)-boundary(j,1))/(boundary(j+4,2)-boundary(j,2))))]; 
% % end
% % 
% % close all;
% % t1=1:size(angle41',2);
% % t2=1:size(angle42',2);
% % t3=1:size(angle43',2);
% % t4=1:size(angle44',2);
% % 
% % display('Std ler');
% % display('Std Angle41');
% % std(angle41)
% % display('Std Angle42');
% % std(angle42)
% % display('Std Angle43');
% % std(angle43)
% % display('Std Angle44');
% % std(angle44)
% % 
% % display('Std lerin Mean');
% % (std(angle41)+std(angle42)+std(angle43)+std(angle44))/4
% % 
% % 
% % 
% % 
% % figure
% % plot(t1,angle41','--rs',...
% %     'LineWidth',1.5,...
% %     'MarkerSize',5,...
% %     'MarkerEdgeColor','g',...
% %     'MarkerFaceColor',[0.3,0.3,0.3])
% % title(['Angles  ' filename]);
% % xlabel('Distance 4 pixels apart');
% % ylabel('Angles');
% % 
% % 
% % figure
% % plot(t2,angle42','--ks',...
% %     'LineWidth',1.5,...
% %     'MarkerSize',5,...
% %     'MarkerEdgeColor','b',...
% %     'MarkerFaceColor',[0.5,0.5,0.5])
% % title(['Absolute Value Angles Averages ' filename]);
% % xlabel('Distance 4 pixels apart');
% % ylabel('Angles');
% % 
% % figure
% % plot(t3,angle43','--bs',...
% %     'LineWidth',1.5,...
% %     'MarkerSize',5,...
% %     'MarkerEdgeColor','m',...
% %     'MarkerFaceColor',[0.5,0.5,0.5])
% % title(['Angles' filename]);
% % xlabel('Distance 4 pixels apart');
% % ylabel('Angles ');
% % 
% % figure
% % plot(t4,angle44','--bs',...
% %     'LineWidth',1.5,...
% %     'MarkerSize',5,...
% %     'MarkerEdgeColor','r',...
% %     'MarkerFaceColor',[0.5,0.5,0.5])
% % title(['Angles' filename]);
% % xlabel('Distance 4 pixels apart');
% % ylabel('Angles ');      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 8 pixels apart

% % % for i=1:8:size(boundary,1)-8
% % %     angle41=[angle41; atand(((boundary(i+8,1)-boundary(i,1))/(boundary(i+8,2)-boundary(i,2))))]; 
% % % end
% % % 
% % % for j=2:8:size(boundary,1)-8
% % %     angle42=[angle42; atand(((boundary(j+8,1)-boundary(j,1))/(boundary(j+8,2)-boundary(j,2))))]; 
% % % end
% % % 
% % % for j=3:8:size(boundary,1)-8
% % %     angle43=[angle43; atand(((boundary(j+8,1)-boundary(j,1))/(boundary(j+8,2)-boundary(j,2))))]; 
% % % end
% % % 
% % % for j=4:8:size(boundary,1)-8
% % %     angle44=[angle44; atand(((boundary(j+8,1)-boundary(j,1))/(boundary(j+8,2)-boundary(j,2))))]; 
% % % end
% % % 
% % % close all;
% % % t1=1:size(angle41',2);
% % % t2=1:size(angle42',2);
% % % t3=1:size(angle43',2);
% % % t4=1:size(angle44',2);
% % % 
% % % display('Std ler');
% % % display('Std Angle81');
% % % std(angle41)
% % % display('Std Angle82');
% % % std(angle42)
% % % display('Std Angle83');
% % % std(angle43)
% % % display('Std Angle84');
% % % std(angle44)
% % % 
% % % display('Std lerin Mean');
% % % (std(angle41)+std(angle42)+std(angle43)+std(angle44))/4
% % % 
% % % 
% % % figure
% % % plot(t1,angle41','--rs',...
% % %     'LineWidth',1.5,...
% % %     'MarkerSize',5,...
% % %     'MarkerEdgeColor','g',...
% % %     'MarkerFaceColor',[0.3,0.3,0.3])
% % % title(['Angles  ' filename]);
% % % xlabel('Distance 4 pixels apart');
% % % ylabel('Angles');
% % % 
% % % 
% % % figure
% % % plot(t2,angle42','--ks',...
% % %     'LineWidth',1.5,...
% % %     'MarkerSize',5,...
% % %     'MarkerEdgeColor','b',...
% % %     'MarkerFaceColor',[0.5,0.5,0.5])
% % % title(['Absolute Value Angles Averages ' filename]);
% % % xlabel('Distance 4 pixels apart');
% % % ylabel('Angles');
% % % 
% % % figure
% % % plot(t3,angle43','--bs',...
% % %     'LineWidth',1.5,...
% % %     'MarkerSize',5,...
% % %     'MarkerEdgeColor','m',...
% % %     'MarkerFaceColor',[0.5,0.5,0.5])
% % % title(['Angles' filename]);
% % % xlabel('Distance 4 pixels apart');
% % % ylabel('Angles ');
% % % 
% % % figure
% % % plot(t4,angle44','--bs',...
% % %     'LineWidth',1.5,...
% % %     'MarkerSize',5,...
% % %     'MarkerEdgeColor','r',...
% % %     'MarkerFaceColor',[0.5,0.5,0.5])
% % % title(['Angles' filename]);
% % % xlabel('Distance 4 pixels apart');
% % % ylabel('Angles ');      




