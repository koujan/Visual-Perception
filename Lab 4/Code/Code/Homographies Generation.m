% Lab 4 - SIFT
% Gourab Ghosh Roy and Mohammad Rami Koujan

clc
clear all
close all
warning off Images:initSize:adjustingMag

I = imread('Image_base_050.jpg');
figure();imshow(I)
Iorig_l = I(2083:3582,2070:3569,:);
figure();imshow(Iorig_l)
Iorig = I(2582:3081,2443:3192,:);
figure();imshow(Iorig)

%% sequence 1

Iorig_w=I(2582:3081,2242:3391,:);
Iorig_w2=I(2381:3280,2443:3192,:);
Sequence1Homographies=struct('H',[]);
if(exist('Sequence1Homographies','file')~=7)
    mkdir('Sequence1Homographies');
end
imwrite(rgb2gray(Iorig),['Sequence1Homographies\Image_00a.png']);
k=50;
X=[1,1;size(Iorig_w,2),1;1,size(Iorig_w,1);size(Iorig_w,2),size(Iorig_w,1)]; % corners of the horizontally stretched image
X2=[1,1;size(Iorig_w2,2),1;1,size(Iorig_w2,1);size(Iorig_w2,2),size(Iorig_w2,1)]; % corners of the vertically stretched image
for i=1:4
    for j=1:4
        U(1:4,:)=[X(1,1)+j*k,X(1,2);X(2,1)-j*k,X(2,2);X(3,1),X(3,2);X(4,1),X(4,2)]; % up
        U(5:8,:)=[X(1,1),X(1,2);X(2,1),X(2,2);X(3,1)+j*k,X(3,2);X(4,1)-j*k,X(4,2)]; % down
        U(9:12,:)=[X2(1,1),X2(1,2);X2(2,1),X2(2,2)+j*k;X2(3,1),X2(3,2);X2(4,1),X2(4,2)-j*k]; % right
        U(13:16,:)=[X2(1,1),X2(1,2)+j*k;X2(2,1),X2(2,2);X2(3,1),X2(3,2)-j*k;X2(4,1),X2(4,2)]; % left
        if(i<3)
             T=fitgeotrans(U((i-1)*4+1:i*4,:),X,'projective');
             Itransform=imwarp(Iorig_w,T);
             Sequence1Homographies(j+(i-1)*4).H= T.T';%T.tdata.T';
             [m,n,~] = size(Itransform);
             Itransform = Itransform(ceil(m/2)-250:ceil(m/2)+249,floor(n/2)-373:floor(n/2)+376,:);
        else
             T=fitgeotrans(U((i-1)*4+1:i*4,:),X2,'projective');
             Itransform=imwarp(Iorig_w2,T);
             Sequence1Homographies(j+(i-1)*4).H=T.T';%T.tdata.T';
             [m,n,~] = size(Itransform);
             Itransform = Itransform(floor(m/2)-248:floor(m/2)+251,floor(n/2)-374:floor(n/2)+375,:);
        end
         imwrite(rgb2gray(Itransform),['Sequence1Homographies\Image_0' num2str(j+(i-1)*4) 'a.png' ]);
         imwrite(rgb2gray(Itransform+uint8(3*randn(size(Itransform)))),['Sequence1Homographies\Image_0' num2str(j+(i-1)*4) 'b.png' ]);
         imwrite(rgb2gray(Itransform+uint8(6*randn(size(Itransform)))),['Sequence1Homographies\Image_0' num2str(j+(i-1)*4) 'c.png' ]);
         imwrite(rgb2gray(Itransform+uint8(18*randn(size(Itransform)))),['Sequence1Homographies\Image_0' num2str(j+(i-1)*4) 'd.png' ]);
         figure,imshow(Itransform);
        
    end
end
save('Sequence1Homographies\Sequence1Homographies.mat','Sequence1Homographies');

%% sequence 2

zoom=1.1:.05:1.5;
Sequence2Homographies=struct('H',[]);
if(exist('Sequence2Homographies','file')~=7)
    mkdir('Sequence2Homographies');
end
imwrite(rgb2gray(Iorig),['Sequence2Homographies\Image_00a.png']);
for i=1:length(zoom)
    Sequence2Homographies(i).H=[1, 0, size(Iorig,2)/2;0 ,-1, size(Iorig,1)/2 ;0 ,0 ,1 ]* [zoom(i) 0 0;0 zoom(i) 0;0,0,1]* [1, 0 , -size(Iorig,2)/2 ; 0 , -1 ,size(Iorig,1)/2; 0, 0, 1];%[zoom(i) 0 -boarder(1,1);0 zoom(i) -boarder(2,1);0,0,1];
    T = maketform('affine',[zoom(i) 0 0; 0 zoom(i) 0; 0 0 1]);
    Itransform = imtransform(Iorig,T);
    [m,n,~] = size(Itransform);
    Itransform = Itransform(round(m/2)-249:round(m/2)+250,round(n/2)-374:round(n/2)+375,:);
    figure();imshow(Itransform);
    imwrite(rgb2gray(Itransform),['Sequence2Homographies\Image_0' num2str(i) 'a.png' ]);
    imwrite(rgb2gray(Itransform+uint8(3*randn(size(Itransform)))),['Sequence2Homographies\Image_0' num2str(i) 'b.png' ]);
    imwrite(rgb2gray(Itransform+uint8(6*randn(size(Itransform)))),['Sequence2Homographies\Image_0' num2str(i) 'c.png' ]);
    imwrite(rgb2gray(Itransform+uint8(18*randn(size(Itransform)))),['Sequence2Homographies\Image_0' num2str(i) 'd.png' ]);
end
save('Sequence2Homographies\Sequence2Homographies.mat','Sequence2Homographies');

%% sequence 3

theta=(-45:5:-5)*pi/180;
theta=[theta, (5:5:45) * pi/180];
Sequence3Homographies=struct('H',[]);
if(exist('Sequence3Homographies','file')~=7)
    mkdir('Sequence3Homographies');
end
imwrite(rgb2gray(Iorig),['Sequence3Homographies\Image_00a.png']);
for i=1:length(theta)
    Sequence3Homographies(i).H=[1, 0, size(Iorig,2)/2;0 ,-1, size(Iorig,1)/2 ;0 ,0 ,1 ]*[cos(theta(i)), -sin(theta(i)), 0; sin(theta(i)) ,cos(theta(i)), 0;0 ,0 ,1] * [1, 0 , -size(Iorig,2)/2 ; 0 , -1 ,size(Iorig,1)/2; 0, 0, 1];
    T = maketform('affine',[cos(theta(i)) -sin(theta(i)) 0; sin(theta(i)) cos(theta(i)) 0; 0 0 1]);
    Itransform =imwarp(Iorig_l,T);
    [m,n,~] = size(Itransform);
    Itransform = Itransform(floor(m/2)-249:floor(m/2)+250,floor(n/2)-374:floor(n/2)+375,:);
    figure();imshow(Itransform);
    imwrite(Itransform,['Sequence3Homographies\Image_0' num2str(i) 'a.png' ]);
    imwrite(Itransform+uint8(3*randn(size(Itransform))),['Sequence3Homographies\Image_0' num2str(i) 'b.png' ]);
    imwrite(Itransform+uint8(6*randn(size(Itransform))),['Sequence3Homographies\Image_0' num2str(i) 'c.png' ]);
    imwrite(Itransform+uint8(18*randn(size(Itransform))),['Sequence3Homographies\Image_0' num2str(i) 'd.png' ]);
    
end
save('Sequence3Homographies\Sequence3Homographies.mat','Sequence3Homographies');

