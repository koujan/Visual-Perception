clear all;clc;
%%%%%%%% part1: %%%%%%%%%%%%%%

% step 1
au=557.0943; av=712.9824; u0=326.3819; v0=298.6679;
f=80;
Tx=100; Ty=0; Tz=1500;
Phix=0.8*pi/2; Phiy=-1.8*pi/2; Phix1=pi/5; %Euler_XYX1 & image size:640*480

% step 2
Rotx=[1 0 0;0 cos(Phix) -sin(Phix);0 sin(Phix) cos(Phix)];
Roty=[cos(Phiy) 0 sin(Phiy) ;0 1 0;-sin(Phiy) 0 cos(Phiy)];
Rotx1=[1 0 0;0 cos(Phix1) -sin(Phix1);0 sin(Phix1) cos(Phix1)];
T=[Tx;Ty;Tz];
R=Rotx*Roty*Rotx1;
ext=[R(1,:) T(1);R(2,:) T(2);R(3,:) T(3);0 0 0 1];
Int=[au 0 u0 0; 0 av v0 0;0 0 1 0];

% step3
Num_points=6; %Num_points=10; Num_points=50; % Number of 3D points
points=randi([-480,480],3,Num_points);

% step4
points=[points;ones(1,Num_points)];
camera_mat=Int*ext;
proj=camera_mat*points;

% step 5
for i=1:size(proj,2)
proj(1:2,i)=proj(1:2,i)/proj(3,i); % normalization by s
end
scatter(proj(1,:),proj(2,:));
close();
% step 6
Q=zeros(size(points,2)*2,11); % every point gives two equations
B=zeros(size(points,2)*2,1);
for i=1:size(proj,2)
Q(2*i-1,:)= [points(1,i),points(2,i), points(3,i),1,0,0,0,0,-proj(1,i)*points(1,i), -proj(1,i)*points(2,i),-proj(1,i)*points(3,i)];    
Q(2*i,:)=[0,0,0,0,points(1,i),points(2,i), points(3,i),1,-proj(2,i)*points(1,i),-proj(2,i)*points(2,i),-proj(2,i)*points(3,i)];
B(2*i-1)=proj(1,i);
B(2*i)=proj(2,i);
end
Hall_mat=Q\B;  % vector of size 11*1
Hall_mat=[Hall_mat; 1];
Hall_mat=[Hall_mat(1:4)';Hall_mat(5:8)';Hall_mat(9:12)'];

% step 7
camera_mat=camera_mat./camera_mat(3,4);
diff=Hall_mat-camera_mat;
ave=sum(diff(:))/12;
SD=sqrt(var(diff(:)));

% step 8
Noise=0.5*randn(size(proj));  % to adapt to the different noise ranges, please change the scaling factor (0.5) accordingly.
proj_n=proj+Noise;
Q_n=zeros(size(points,2)*2,11);
B_n=zeros(size(points,2)*2,1);
for i=1:size(proj,2)
Q_n(2*i-1,:)= [points(1,i),points(2,i), points(3,i),1,0,0,0,0,-proj_n(1,i)*points(1,i), -proj_n(1,i)*points(2,i),-proj_n(1,i)*points(3,i)];    
Q_n(2*i,:)=[0,0,0,0,points(1,i),points(2,i), points(3,i),1,-proj_n(2,i)*points(1,i),-proj_n(2,i)*points(2,i),-proj_n(2,i)*points(3,i)];
B_n(2*i-1)=proj_n(1,i);
B_n(2*i)=proj_n(2,i);
end
Hall_mat_n=Q_n\B_n;
Hall_mat_n=[Hall_mat_n; 1];
Hall_mat_n=[Hall_mat_n(1:4)';Hall_mat_n(5:8)';Hall_mat_n(9:12)'];
diff_n=Hall_mat_n-Hall_mat;
ave_n=sum(diff_n(:))/12;
SD_n=sqrt(var(diff_n(:)));
proj_n_r=Hall_mat_n*points;
for i=1:size(proj,2)
proj_n_r(1:2,i)=proj_n_r(1:2,i)/proj_n_r(3,i); % normalization by s
end
dis_norm=zeros(1,Num_points);
for i=1:size(proj,2)
dis_norm(i)=norm([(proj(1,i)-proj_n_r(1,i)),(proj(2,i)-proj_n_r(2,i))]);
end
ave_n2=mean(dis_norm);
SD_n2=sqrt(var(dis_norm));

%% part2

% step 10
Q_faug=zeros(size(points,2)*2,11);
B_faug=zeros(size(points,2)*2,1);
for i=1:size(proj,2)
Q_faug(2*i-1,:)= [points(1,i),points(2,i), points(3,i),-proj(1,i)*points(1,i),-proj(1,i)*points(2,i),-proj(1,i)*points(3,i),0,0,0,1,0];
Q_faug(2*i,:)=[0,0,0,-proj(2,i)*points(1,i),-proj(2,i)*points(2,i),-proj(2,i)*points(3,i),points(1,i),points(2,i), points(3,i),0,1];
B_faug(2*i-1)=proj(1,i);
B_faug(2*i)=proj(2,i);
end
X=Q_faug\B_faug;
T1=X(1:3);
T2=X(4:6);
T3=X(7:9);
C1=X(10);
C2=X(11);
u01=dot(T1,T2)/(norm(T2)^2);
v01=dot(T2,T3)/(norm(T2)^2);
au1=norm(cross(T1,T2))/(norm(T2)^2);
av1=norm(cross(T2,T3))/(norm(T2)^2);
r1=( norm(T2)/norm(cross(T1,T2)) )*( T1-( ( dot(T1,T2)/(norm(T2)^2) )*T2 ) );
r2=( norm(T2)/norm(cross(T2,T3)) )*( T3-( ( dot(T2,T3)/(norm(T2)^2) )*T2 ) );
r3=T2/norm(T2);
tx=( norm(T2)/norm(cross(T1,T2)) )*( C1-( dot(T1,T2)/(norm(T2)^2) ));
ty=( norm(T2)/norm(cross(T2,T3)) )*( C2-( dot(T2,T3)/(norm(T2)^2) ));
tz=1/norm(T2);
faug_mat=[T1'*tz C1*tz;T3'*tz C2*tz;T2'*tz tz];
R_diff=R-[r1';r2';r3'];
T_diff=[Tx;Ty;Tz]-[tx;ty;tz];
u0_diff=u01-u0;
v0_diff=v01-v0;
au_diff=au1-au;
av_diff=av1-av;

% step 11
proj_nn1=proj_n;%proj_n is defined in step 8
Q_n1=zeros(size(points,2)*2,11);
B_n1=zeros(size(points,2)*2,1);
for i=1:size(proj,2)
Q_n1(2*i-1,:)= [points(1,i),points(2,i), points(3,i),-proj_nn1(1,i)*points(1,i),-proj_nn1(1,i)*points(2,i),-proj_nn1(1,i)*points(3,i),0,0,0,1,0];
Q_n1(2*i,:)=[0,0,0,-proj_nn1(2,i)*points(1,i),-proj_nn1(2,i)*points(2,i),-proj_nn1(2,i)*points(3,i),points(1,i),points(2,i), points(3,i),0,1];
B_n1(2*i-1)=proj_nn1(1,i);
B_n1(2*i)=proj_nn1(2,i);
end
X_n1=Q_n1\B_n1;
T1=X_n1(1:3);
T2=X_n1(4:6);
T3=X_n1(7:9);
C1=X_n1(10);
C2=X_n1(11);
tz=1/norm(T2);
faug_mat_n1=[T1' C1;T3' C2;T2' 1]; % normalized by tz to compare it with Hall's
%faug_mat_n1=[T1'*tz C1*tz;T3'*tz C2*tz;T2'*tz tz];
proj_nn1_r=faug_mat_n1*points;
for i=1:size(proj_nn1_r,2)
proj_nn1_r(1:2,i)=proj_nn1_r(1:2,i)/proj_nn1_r(3,i); % normalization by s
end
dis_norm_nn1=zeros(1,Num_points);
for i=1:size(proj,2)
dis_norm_nn1(i)=norm([(proj(1,i)-proj_nn1_r(1,i)),(proj(2,i)-proj_nn1_r(2,i))]);
end
ave_nn1=mean(dis_norm_nn1); % mean of distances 
SD_nn1=sqrt(var(dis_norm_nn1)) ;% standard deviation of distances

% same accuracy for Hall and Faugeras points;


%% Part 3
W=3000;
% world x axis
scatter3(points(1,:),points(2,:),points(3,:),'fill');hold on;
line([W 0],[0 0],[0 0]);
line([W W-100],[0 50],[0 0]); % first part of the arrow
line([W W-100],[0 -50],[0 0]);% second part of the arrow
text(W,0,0,'X axis');
% world y axis
line([0 0],[0 W],[0 0]);
line([0 50],[W W-100],[0 0]);% first part of the arrow
line([0 -50],[W W-100],[0 0]);% second part of the arrow
text(0,W,0,'Y axis');
% world z axis
line([0 0],[0 0],[0 W]);
line([0 0],[0 -50],[W W-100]);% first part of the arrow
line([0 0],[0 +50],[W W-100]);% second part of the arrow
text(0,0,W,'Z axis');
text(0,0,0,'Ow');

% calculating the inverse transformation matrix
inv_ext=[R',-R'*T];
inv_ext=[inv_ext;0 0 0 1];
% inv_ext*ext: should be equal to identity

% computing the origin and standard unit vectors of the camersa coordinate system
Oc=inv_ext*[0;0;0;1];
Xc=inv_ext*[W;0;0;1];
Yc=inv_ext*[0;W;0;1];
Zc=inv_ext*[0;0;W;1];

% camera x axis
line([Oc(1) Xc(1)],[Oc(2) Xc(2)],[Oc(3) Xc(3)]);hold on;
text(Xc(1),Xc(2),Xc(3),'X axis');

% camera y axis
line([Oc(1) Yc(1)],[Oc(2) Yc(2)],[Oc(3) Yc(3)]);
text(Yc(1),Yc(2),Yc(3),'Y axis');

% camera z axis
line([Oc(1) Zc(1)],[Oc(2) Zc(2)],[Oc(3) Zc(3)]);
text(Zc(1),Zc(2),Zc(3),'Z axis');
text(Oc(1),Oc(2),Oc(3),'Oc');

% checking that camera unit vectors are perpendicular
%dot([Xc(1:3)-Oc(1:3)],[Yc(1:3)-Oc(1:3)])
%dot([Zc(1:3)-Oc(1:3)],[Yc(1:3)-Oc(1:3)])
%dot([Zc(1:3)-Oc(1:3)],[Xc(1:3)-Oc(1:3)])

% focal point
foc=inv_ext*[0;0;f;1];
scatter3(foc(1),foc(2),foc(3),'k','fill');


% image plane
Ix=u0*f/au;
Iy=v0*f/av; 
im_cor1=[-Ix,-Iy,f];im_cor1=inv_ext*[im_cor1,1]';
im_cor2=[Ix,-Iy,f];im_cor2=inv_ext*[im_cor2,1]';
im_cor3=[Ix,Iy,f];im_cor3=inv_ext*[im_cor3,1]';
im_cor4=[-Ix,Iy,f];im_cor4=inv_ext*[im_cor4,1]';
line([im_cor1(1) im_cor2(1)],[im_cor1(2) im_cor2(2)],[im_cor1(3) im_cor2(3)],'Color','g');
line([im_cor2(1) im_cor3(1)],[im_cor2(2) im_cor3(2)],[im_cor2(3) im_cor3(3)],'Color','g');
line([im_cor3(1) im_cor4(1)],[im_cor3(2) im_cor4(2)],[im_cor3(3) im_cor4(3)],'Color','g');
line([im_cor4(1) im_cor1(1)],[im_cor4(2) im_cor1(2)],[im_cor4(3) im_cor1(3)],'Color','g');
%text(im_cor4(1),im_cor4(2),im_cor4(3),'Image Plane');


% 2d points
points_2d=ext*points; % refrences with respect to camera
points_2d(1,:)=f*points_2d(1,:)./points_2d(3,:);
points_2d(2,:)=f*points_2d(2,:)./points_2d(3,:);
points_2d(3,:)=f;
points_2d=inv_ext*points_2d;
scatter3(points_2d(1,:),points_2d(2,:),points_2d(3,:));

%lines between the 3d points and their projections on the image plane
for i=1:size(points,2)
   line([Oc(1,1) points(1,i)],[Oc(2,1) points(2,i)],[Oc(3,1) points(3,i)]); 
end

% checking step to make sure that each 3d point with its projection,2d point , and the origin of the camera coordinate system
% are collinear. The result of the cross prodcut should be theoritically zero:
for i=1:size(points,2)
    cross(points_2d(1:3,i)-Oc(1:3),points(1:3,i)-Oc(1:3));
end


