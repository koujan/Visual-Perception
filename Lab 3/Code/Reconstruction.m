clear all;clc
%% step 1
au1 = 100; av1 = 120; uo1 = 128; vo1 = 128; % image size: 256*256

%% step 2
au2 = 90; av2 = 110; uo2 = 128; vo2 = 128;
ax = 0.1; by = pi/4; cz = 0.2 ;    % XYZ Euler
tx = -1000; ty = 190; tz = 230;    % image size: 256*256
Rotx=[1 0 0;0 cos(ax) -sin(ax);0 sin(ax) cos(ax)];
Roty=[cos(by) 0 sin(by) ;0 1 0;-sin(by) 0 cos(by)]; 
Rotz=[cos(cz) -sin(cz) 0;sin(cz) cos(cz) 0;0 0 1];
T=[tx,ty,tz];
R=Rotx*Roty*Rotz;

%% stpe 3
p1=[eye(3) [0 0 0]'];
p2=[R T'];
A1=[au1 0 uo1;0 av1 vo1;0 0 1];
A2=[au2 0 uo2;0 av2 vo2;0 0 1];

%% step 4
antiSym=[0 -T(3) T(2); T(3) 0 -T(1); -T(2) T(1) 0];
F=inv(A2')*R'*antiSym*inv(A1);

%% step 5
V(:,1) = [100;-400;2000;1];
V(:,2) = [300;-400;3000;1];
V(:,3) = [500;-400;4000;1];
V(:,4) = [700;-400;2000;1];
V(:,5) = [900;-400;3000;1];
V(:,6) = [100;-50;4000;1];
V(:,7) = [300;-50;2000;1];
V(:,8) = [500;-50;3000;1];
V(:,9) = [700;-50;4000;1];
V(:,10) = [900;-50;2000;1];
V(:,11) = [100;50;3000;1];
V(:,12) = [300;50;4000;1];
V(:,13) = [500;50;2000;1];
V(:,14) = [700;50;3000;1];
V(:,15) = [900;50;4000;1];
V(:,16) = [100;400;2000;1];
V(:,17) = [300;400;3000;1];
V(:,18) = [500;400;4000;1];
V(:,19) = [700;400;2000;1];
V(:,20) = [900;400;3000;1];
% xx=randi([100,900],1,30);
% yy=randi([-400,400],1,30);
% zz=randi([2000,4000],1,30);
% V(:,21:50)=[xx;yy;zz;ones(1,30)];
%% step 6
Int1=[A1,[0 0 0]'];
Int2=[A2,[0 0 0]'];
ext1=[p1;[0 0 0 1]];
ext2=[R' -R'*T';0 0 0 1];
points1=Int1*ext1*V;
points2=Int2*ext2*V;
for i=1:size(V,2)
    points1(1:3,i)=points1(1:3,i)/points1(3,i);
    points2(1:3,i)=points2(1:3,i)/points2(3,i);
end

%% step 7
%fig1=figure;scatter(points1(1,:),points1(2,:));
fig2=figure;scatter(points2(1,:),points2(2,:));hold on;

%% step 8
u=zeros(size(V,2),8);
for i=1:size(V,2)
    u(i,:)=[points1(1,i)*points2(1,i),points1(2,i)*points2(1,i), points2(1,i),points1(1,i)*points2(2,i),points1(2,i)*points2(2,i),points2(2,i),points1(1,i),points1(2,i)];
end
F_least=-u\ones(size(V,2),1);
F_least=[F_least ;1];

%% step 9
F_least=[F_least(1:3)';F_least(4:6)';F_least(7:9)'];
F=F(:,:)/F(3,3);
F_diff=F-F_least;

%% step 10

%%%%%% epipole geometry for the second image
L_prime=F_least*points1;
x1=-500;
x2=500;
m=zeros(size(V,2),1);
d=zeros(size(V,2),1);
for i=1:size(V,2)
m(i)=-L_prime(1,i)/L_prime(2,i);
d(i)=-L_prime(3,i)/L_prime(2,i);
y1=m(i)*x1+d(i);
y2=m(i)*x2+d(i);

plot([x1 x2],[y1 y2]);hold on;
end

% computing the epipole by projection
e_prime=Int2*ext2*[0;0;0;1];
e_prime=e_prime(1:2,1)/e_prime(3,1);
scatter(e_prime(1),e_prime(2),'fill');
text(e_prime(1),e_prime(2),'  Epipole\_prime');
title('Second image plane');

% computing the epipole from the crossing of the epipolar lines
mat=[-m,ones(size(V,2),1)];
e_prime_cross=mat\d;
e_diff_cross=e_prime-e_prime_cross;

% computing the epipole from the Fundamental matrix
[~,~,v_mat] = svd(F');
e_prime_Fund=v_mat(1:2,end)/v_mat(end,end);
e_diff_Fund=e_prime-e_prime_Fund;


%%%%%%%%%% eipipole geometry for first image plane
L=F_least'*points2;
figure;scatter(points1(1,:),points1(2,:));hold on;
m2=zeros(size(V,2),1);
d2=zeros(size(V,2),1);
for i=1:size(V,2)
m2(i)=-L(1,i)/L(2,i);
d2(i)=-L(3,i)/L(2,i);
y1=m2(i)*x1+d2(i);
y2=m2(i)*x2+d2(i);
plot([x1 x2],[y1 y2]);
end

% computing the epipole by projection
e=Int1*ext1*[T,1]';
e=e(1:2)/e(3);
scatter(e(1),e(2),'fill');
text(e(1),e(2),'  Epipole');
title('First image plane');

% computing the epipole from the crossing of the epipolar lines
mat2=[-m2,ones(size(V,2),1)];
e_cross=mat2\d2;
e_diff2_cross=e-e_cross;

% computing the epipole from the Fundamental matrix
[~,~,v_mat] = svd(F);
e_Fund=v_mat(1:2,end)/v_mat(end,end);
e_diff2_Fund=e-e_Fund;

%% Step 11
sigma=0.5;
Noise=sigma*randn(2,size(points1,2));
Noise2=sigma*randn(2,size(points2,2));
points1_n=points1(1:2,:)+Noise;points1_n(3,:)=1;
points2_n=points2(1:2,:)+Noise2;points2_n(3,:)=1;

%% step 12

u_n=zeros(size(V,2),8);
for i=1:size(V,2)
    u_n(i,:)=[points1_n(1,i)*points2_n(1,i),points1_n(2,i)*points2_n(1,i), points2_n(1,i),points1_n(1,i)*points2_n(2,i),points1_n(2,i)*points2_n(2,i),points2_n(2,i),points1_n(1,i),points1_n(2,i)];
end
F_n=-u_n\ones(size(V,2),1);
F_n=[F_n ;1];

%step 12-9
F_n=[F_n(1:3)';F_n(4:6)';F_n(7:9)'];
F_diff_n=F_n-F;

%step 12-10

%%%%%% epipole geometry for the second image
L_prime_n=F_n*points1_n;
figure;scatter(points2_n(1,:),points2_n(2,:));hold on;
m_n=zeros(size(V,2),1);
d_n=zeros(size(V,2),1);
mean_dis=zeros(size(V,2),1);
for i=1:size(V,2)
m_n(i)=-L_prime_n(1,i)/L_prime_n(2,i);
d_n(i)=-L_prime_n(3,i)/L_prime_n(2,i);
y1=m_n(i)*x1+d_n(i);
y2=m_n(i)*x2+d_n(i);
mean_dis(i)=abs(L_prime_n(1,i)*points2_n(1,i)+L_prime_n(2,i)*points2_n(2,i)+L_prime_n(3,i))/sqrt(L_prime_n(1,i)^2+L_prime_n(2,i)^2);
plot([x1 x2],[y1 y2]);hold on;
end
ave=mean(mean_dis);
s_dev=std(mean_dis);

% computing the epipole by projection    
        %  the same as without noise  %
scatter(e_prime(1),e_prime(2),'fill');
text(e_prime(1),e_prime(2),'  Epipole\_prime');
title('Second image plane (noisy)');
 
% computing the epipole from the crossing of the epipolar lines
mat_n=[-m_n,ones(size(V,2),1)];
e_prime_cross_n=mat_n\d_n;
e_diff_cross_n=e_prime(1:2)-e_prime_cross_n;
scatter(e_prime_cross_n(1),e_prime_cross_n(2),'fill');
text(e_prime_cross_n(1),e_prime_cross_n(2),'  Epipole\_prime\_cross\_noisy');

% computing the epipole from the fundamental matrix
[~,~,v_mat_n] = svd(F_n');
e_prime_Fund_n=v_mat_n(1:2,end)/v_mat_n(end,end);
e_diff_Fund_n=e_prime-e_prime_Fund_n;
scatter(e_prime_Fund_n(1),e_prime_Fund_n(2),'fill');
text(e_prime_Fund_n(1),e_prime_Fund_n(2),'  Epipole\_prime\_F\_noisy');


%%%%%%%%%%% eipipole geometry for first image plane
L_n=F_n'*points2_n;
figure;scatter(points1_n(1,:),points1_n(2,:));hold on;
m2_n=zeros(size(V,2),1);
d2_n=zeros(size(V,2),1);
mean_dis2=zeros(size(V,2),1);
for i=1:size(V,2)
m2_n(i)=-L_n(1,i)/L_n(2,i);
d2_n(i)=-L_n(3,i)/L_n(2,i);
y1=m2_n(i)*x1+d2_n(i);
y2=m2_n(i)*x2+d2_n(i);
mean_dis2(i)=abs(L_n(1,i)*points1_n(1,i)+L_n(2,i)*points1_n(2,i)+L_n(3,i))/sqrt(L_n(1,i)^2+L_n(2,i)^2);
%mean_dis2(i)=abs(det([x2-x1,points1_n(1,i)-x1;y2-y1,points1_n(2,i)-y1]))/norm([x2-x1,y2-y1]);
plot([x1 x2],[y1 y2]);
end
ave2=mean(mean_dis2);
s_dev2=std(mean_dis2);
ave_LS=(ave+ave2)/2
s_dev_LS=(s_dev+s_dev2)/2

% computing the epipole by projection
      
        %  the same as without noise  %
title('First image plane (noisy)');
scatter(e(1),e(2),'fill');
text(e(1),e(2),'  Epipole');

% computing the epipole from the crossing of the epipolar lines
mat2_n=[-m2_n,ones(size(V,2),1)];
e_cross_n=mat2_n\d2_n;
e_diff2_n=e(1:2)-e_cross_n;
scatter(e_cross_n(1),e_cross_n(2),'fill');
text(e_cross_n(1),e_cross_n(2),'  Epipole\_cross\_noisy');

% computing the epipole from the fundamental matrix
[~,~,v_mat_n] = svd(F_n);
e_Fund_n=v_mat_n(1:2,end)/v_mat_n(end,end);
e_diff_Fund_n2=e-e_Fund_n;
scatter(e_Fund_n(1),e_Fund_n(2),'fill');
text(e_Fund_n(1),e_Fund_n(2),'  Epipole\_F\_noisy');




                                      %%%%%%%%    SVD     %%%%%%%%%%
[U,S,v_mat] = svd(F_n);
[~,i]=min([S(1,1),S(2,2),S(3,3)]);
S(i,i)=0;
F_n=U*S*v_mat';


% draw again

%step 12-10

%%%%%% epipole geometry for the second image
L_prime_n=F_n*points1_n;
figure;scatter(points2_n(1,:),points2_n(2,:));hold on;
m_n=zeros(size(V,2),1);
d_n=zeros(size(V,2),1);
for i=1:size(V,2)
m_n(i)=-L_prime_n(1,i)/L_prime_n(2,i);
d_n(i)=-L_prime_n(3,i)/L_prime_n(2,i);
y1=m_n(i)*x1+d_n(i);
y2=m_n(i)*x2+d_n(i);
plot([x1 x2],[y1 y2]);hold on;
end

% computing the epipole by projection    
        %  the same as without noise  %
scatter(e_prime(1),e_prime(2),'fill');
text(e_prime(1),e_prime(2),'  Epipole\_prime');
title('Second image plane (noisy & SVD)');
 
% computing the epipole from the crossing of the epipolar lines
mat_n=[-m_n,ones(size(V,2),1)];
e_prime_cross_n=mat_n\d_n;
e_diff_cross_n=e_prime(1:2)-e_prime_cross_n;
scatter(e_prime_cross_n(1),e_prime_cross_n(2),'fill');
text(e_prime_cross_n(1),e_prime_cross_n(2),'  Epipole\_prime\_cross\_noisy');

% computing the epipole from the fundamental matrix
[~,~,v_mat_n] = svd(F_n');
e_prime_Fund_n=v_mat_n(1:2,end)/v_mat_n(end,end);
e_diff_Fund_n=e_prime-e_prime_Fund_n;
scatter(e_prime_Fund_n(1),e_prime_Fund_n(2),'fill');
text(e_prime_Fund_n(1),e_prime_Fund_n(2),'  Epipole\_prime\_F\_noisy');


%%%%%%%%%%% eipipole geometry for first image plane
L_n=F_n'*points2_n;
figure;scatter(points1_n(1,:),points1_n(2,:));hold on;
m2_n=zeros(size(V,2),1);
d2_n=zeros(size(V,2),1);
for i=1:size(V,2)
m2_n(i)=-L_n(1,i)/L_n(2,i);
d2_n(i)=-L_n(3,i)/L_n(2,i);
y1=m2_n(i)*x1+d2_n(i);
y2=m2_n(i)*x2+d2_n(i);
plot([x1 x2],[y1 y2]);
end
% computing the epipole by projection
      
        %  the same as without noise  %
title('First image plane (noisy & SVD)');
scatter(e(1),e(2),'fill');
text(e(1),e(2),'  Epipole');

% computing the epipole from the crossing of the epipolar lines
mat2_n=[-m2_n,ones(size(V,2),1)];
e_cross_n=mat2_n\d2_n;
e_diff2_n=e(1:2)-e_cross_n;
scatter(e_cross_n(1),e_cross_n(2),'fill');
text(e_cross_n(1),e_cross_n(2),'  Epipole\_cross\_noisy');

% computing the epipole from the fundamental matrix
[~,~,v_mat_n] = svd(F_n);
e_Fund_n=v_mat_n(1:2,end)/v_mat_n(end,end);
e_diff_Fund_n2=e-e_Fund_n;
scatter(e_Fund_n(1),e_Fund_n(2),'fill');
text(e_Fund_n(1),e_Fund_n(2),'  Epipole\_F\_noisy');
%close all;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%% step 14  (part2)
u_SVD=zeros(size(V,2),9);
for i=1:size(V,2)
    u_SVD(i,:)=[points1(1,i)*points2(1,i),points1(2,i)*points2(1,i), points2(1,i),points1(1,i)*points2(2,i),points1(2,i)*points2(2,i),points2(2,i),points1(1,i),points1(2,i),1];
end
[~,S,v_mat]=svd(u_SVD);
[~,i]=min([S(1,1),S(2,2),S(3,3),S(4,4),S(5,5),S(6,6),S(7,7),S(8,8),S(9,9)]);
F_SVD=v_mat(:,i);
F_SVD=F_SVD/F_SVD(end);
F_SVD=[F_SVD(1:3)';F_SVD(4:6)';F_SVD(7:9)'];
F_diff=F_SVD-F_least;


%% step 10 (part 2 repeating steps 10-13)

%%%%%% epipole geometry for the second image
L_prime=F_SVD*points1;
x1=-500;
x2=500;
m=zeros(size(V,2),1);
d=zeros(size(V,2),1);
for i=1:size(V,2)
m(i)=-L_prime(1,i)/L_prime(2,i);
d(i)=-L_prime(3,i)/L_prime(2,i);
y1=m(i)*x1+d(i);
y2=m(i)*x2+d(i);

plot([x1 x2],[y1 y2]);hold on;
end

% computing the epipole by projection
e_prime=Int2*ext2*[0;0;0;1];
e_prime=e_prime(1:2,1)/e_prime(3,1);
figure,scatter(e_prime(1),e_prime(2),'fill');
text(e_prime(1),e_prime(2),'  Epipole\_prime');
title('Second image plane (Part 2)');

% computing the epipole from the crossing of the epipolar lines
mat=[-m,ones(size(V,2),1)];
e_prime_cross=mat\d;
e_diff_cross=e_prime-e_prime_cross;

% computing the epipole from the Fundamental matrix
[~,~,v_mat] = svd(F_SVD');
e_prime_Fund=v_mat(1:2,end)/v_mat(end,end);
e_diff_Fund=e_prime-e_prime_Fund;


%%%%%%%%%% eipipole geometry for first image plane
L=F_SVD'*points2;
figure;scatter(points1(1,:),points1(2,:));hold on;
m2=zeros(size(V,2),1);
d2=zeros(size(V,2),1);
for i=1:size(V,2)
m2(i)=-L(1,i)/L(2,i);
d2(i)=-L(3,i)/L(2,i);
y1=m2(i)*x1+d2(i);
y2=m2(i)*x2+d2(i);
plot([x1 x2],[y1 y2]);
end

% computing the epipole by projection
e=Int1*ext1*[T,1]';
e=e(1:2)/e(3);
scatter(e(1),e(2),'fill');
text(e(1),e(2),'  Epipole');
title('First image plane (Part 2)');

% computing the epipole from the crossing of the epipolar lines
mat2=[-m2,ones(size(V,2),1)];
e_cross=mat2\d2;
e_diff2_cross=e-e_cross;

% computing the epipole from the Fundamental matrix
[~,~,v_mat] = svd(F_SVD);
e_Fund=v_mat(1:2,end)/v_mat(end,end);
e_diff2_Fund=e-e_Fund;


%% Step 11  (part 2 repeating steps 10-13)
% Same as the noise of the LS method defined above in order to compare both results

%% step 12  (part 2 repeating steps 10-13)

u_SVD_n=zeros(size(V,2),9);
for i=1:size(V,2)
    u_SVD_n(i,:)=[points1_n(1,i)*points2_n(1,i),points1_n(2,i)*points2_n(1,i), points2_n(1,i),points1_n(1,i)*points2_n(2,i),points1_n(2,i)*points2_n(2,i),points2_n(2,i),points1_n(1,i),points1_n(2,i),1];
end
[~,S,v_mat]=svd(u_SVD_n);
[~,i]=min([S(1,1),S(2,2),S(3,3),S(4,4),S(5,5),S(6,6),S(7,7),S(8,8),S(9,9)]);
F_SVD_n=v_mat(:,i);
F_SVD_n=F_SVD_n/F_SVD_n(end);
F_SVD_n=[F_SVD_n(1:3)';F_SVD_n(4:6)';F_SVD_n(7:9)'];
F_diff_n=F_SVD_n-F_SVD;

%step 12-10

%%%%%% epipole geometry for the second image
L_prime_n=F_SVD_n*points1_n;
figure;scatter(points2_n(1,:),points2_n(2,:));hold on;
m_n=zeros(size(V,2),1);
d_n=zeros(size(V,2),1);
mean_dis=zeros(size(V,2),1);
for i=1:size(V,2)
m_n(i)=-L_prime_n(1,i)/L_prime_n(2,i);
d_n(i)=-L_prime_n(3,i)/L_prime_n(2,i);
y1=m_n(i)*x1+d_n(i);
y2=m_n(i)*x2+d_n(i);
mean_dis(i)=abs(L_prime_n(1,i)*points2_n(1,i)+L_prime_n(2,i)*points2_n(2,i)+L_prime_n(3,i))/sqrt(L_prime_n(1,i)^2+L_prime_n(2,i)^2);
plot([x1 x2],[y1 y2]);hold on;
end
ave=mean(mean_dis);
s_dev=std(mean_dis);


% computing the epipole by projection    
        %  the same as without noise  %
scatter(e_prime(1),e_prime(2),'fill');
text(e_prime(1),e_prime(2),'  Epipole\_prime');
title('Second image plane (noisy) ,Part2');
 
% computing the epipole from the crossing of the epipolar lines
mat_n=[-m_n,ones(size(V,2),1)];
e_prime_cross_n=mat_n\d_n;
e_diff_cross_n=e_prime(1:2)-e_prime_cross_n;
scatter(e_prime_cross_n(1),e_prime_cross_n(2),'fill');
text(e_prime_cross_n(1),e_prime_cross_n(2),'  Epipole\_prime\_cross\_noisy');

% computing the epipole from the fundamental matrix
[~,~,v_mat_n] = svd(F_SVD_n');
e_prime_Fund_n=v_mat_n(1:2,end)/v_mat_n(end,end);
e_diff_Fund_n=e_prime-e_prime_Fund_n;
scatter(e_prime_Fund_n(1),e_prime_Fund_n(2),'fill');
text(e_prime_Fund_n(1),e_prime_Fund_n(2),'  Epipole\_prime\_F\_noisy');


%%%%%%%%%%% eipipole geometry for first image plane
L_n=F_SVD_n'*points2_n;
figure;scatter(points1_n(1,:),points1_n(2,:));hold on;
m2_n=zeros(size(V,2),1);
d2_n=zeros(size(V,2),1);
mean_dis2=zeros(size(V,2),1);
for i=1:size(V,2)
m2_n(i)=-L_n(1,i)/L_n(2,i);
d2_n(i)=-L_n(3,i)/L_n(2,i);
y1=m2_n(i)*x1+d2_n(i);
y2=m2_n(i)*x2+d2_n(i);
mean_dis2(i)=abs(L_n(1,i)*points1_n(1,i)+L_n(2,i)*points1_n(2,i)+L_n(3,i))/sqrt(L_n(1,i)^2+L_n(2,i)^2);
plot([x1 x2],[y1 y2]);
end
ave2=mean(mean_dis2);
s_dev2=std(mean_dis2);
ave_SVD=(ave+ave2)/2
s_dev_SVD=(s_dev+s_dev2)/2

% computing the epipole by projection
      
        %  the same as without noise  %
title('First image plane (noisy), Part 2');
scatter(e(1),e(2),'fill');
text(e(1),e(2),'  Epipole');

% computing the epipole from the crossing of the epipolar lines
mat2_n=[-m2_n,ones(size(V,2),1)];
e_cross_n=mat2_n\d2_n;
e_diff2_n=e(1:2)-e_cross_n;
scatter(e_cross_n(1),e_cross_n(2),'fill');
text(e_cross_n(1),e_cross_n(2),'  Epipole\_cross\_noisy');

% computing the epipole from the fundamental matrix
[~,~,v_mat_n] = svd(F_SVD_n);
e_Fund_n=v_mat_n(1:2,end)/v_mat_n(end,end);
e_diff_Fund_n2=e-e_Fund_n;
scatter(e_Fund_n(1),e_Fund_n(2),'fill');
text(e_Fund_n(1),e_Fund_n(2),'  Epipole\_F\_noisy');



                            %%%%%%%%    SVD      %%%%%%%%%%
[U,S,v_mat] = svd(F_SVD_n);
[~,i]=min([S(1,1),S(2,2),S(3,3)]);
S(i,i)=0;
F_SVD_n=U*S*v_mat';


% draw again

%step 12-10

%%%%%% epipole geometry for the second image
L_prime_n=F_SVD_n*points1_n;
figure;scatter(points2_n(1,:),points2_n(2,:));hold on;
m_n=zeros(size(V,2),1);
d_n=zeros(size(V,2),1);
for i=1:size(V,2)
m_n(i)=-L_prime_n(1,i)/L_prime_n(2,i);
d_n(i)=-L_prime_n(3,i)/L_prime_n(2,i);
y1=m_n(i)*(x1-500)+d_n(i);
y2=m_n(i)*x2+d_n(i);
plot([x1-500 x2],[y1 y2]);hold on;
end

% computing the epipole by projection    
        %  the same as without noise  %
scatter(e_prime(1),e_prime(2),'fill');
text(e_prime(1),e_prime(2),'  Epipole\_prime');
title('Second image plane (noisy & SVD), Part 2');
 
% computing the epipole from the crossing of the epipolar lines
mat_n=[-m_n,ones(size(V,2),1)];
e_prime_cross_n=mat_n\d_n;
e_diff_cross_n=e_prime(1:2)-e_prime_cross_n;
scatter(e_prime_cross_n(1),e_prime_cross_n(2),'fill');
text(e_prime_cross_n(1),e_prime_cross_n(2),'  Epipole\_prime\_cross\_noisy');

% computing the epipole from the fundamental matrix
[~,~,v_mat_n] = svd(F_SVD_n');
e_prime_Fund_n=v_mat_n(1:2,end)/v_mat_n(end,end);
e_diff_Fund_n=e_prime-e_prime_Fund_n;
scatter(e_prime_Fund_n(1),e_prime_Fund_n(2),'fill');
text(e_prime_Fund_n(1),e_prime_Fund_n(2),'  Epipole\_prime\_F\_noisy');


%%%%%%%%%%% eipipole geometry for first image plane
L_n=F_SVD_n'*points2_n;
figure;scatter(points1_n(1,:),points1_n(2,:));hold on;
m2_n=zeros(size(V,2),1);
d2_n=zeros(size(V,2),1);
for i=1:size(V,2)
m2_n(i)=-L_n(1,i)/L_n(2,i);
d2_n(i)=-L_n(3,i)/L_n(2,i);
y1=m2_n(i)*x1+d2_n(i);
y2=m2_n(i)*x2+d2_n(i);
plot([x1 x2],[y1 y2]);
end
% computing the epipole by projection
      
        %  the same as without noise  %
title('First image plane (noisy & SVD), Part 2');
scatter(e(1),e(2),'fill');
text(e(1),e(2),'  Epipole');

% computing the epipole from the crossing of the epipolar lines
mat2_n=[-m2_n,ones(size(V,2),1)];
e_cross_n=mat2_n\d2_n;
e_diff2_n=e(1:2)-e_cross_n;
scatter(e_cross_n(1),e_cross_n(2),'fill');
text(e_cross_n(1),e_cross_n(2),'  Epipole\_cross\_noisy');

% computing the epipole from the fundamental matrix
[~,~,v_mat_n] = svd(F_SVD_n);
e_Fund_n=v_mat_n(1:2,end)/v_mat_n(end,end);
e_diff_Fund_n2=e-e_Fund_n;
scatter(e_Fund_n(1),e_Fund_n(2),'fill');
text(e_Fund_n(1),e_Fund_n(2),'  Epipole\_F\_noisy');
%close all;
'difference between both methods (LS-SVD)'
ave_LS-ave_SVD
s_dev_LS-s_dev_SVD


%% Optional part

W=3000;
% world x axis
figure,scatter3(V(1,:),V(2,:),V(3,:),'fill');hold on;
line([W 0],[0 0],[0 0]);
line([W W-100],[0 50],[0 0]);  % first part of the arrow
line([W W-100],[0 -50],[0 0]); % second part of the arrow
 text(W,0,0,'X axis');
% % world y axis
line([0 0],[0 W],[0 0]);
line([0 50],[W W-100],[0 0]);% first part of the arrow
line([0 -50],[W W-100],[0 0]);% second part of the arrow
text(0,W,0,'Y axis');
% % world z axis
line([0 0],[0 0],[0 W]);
line([0 0],[0 -50],[W W-100]);% first part of the arrow
line([0 0],[0 +50],[W W-100]);% second part of the arrow
text(0,0,W,'Z axis');
text(0,0,0,'Oc1');

% % focal point
f=40;
foc=[0;0;f];
scatter3(foc(1),foc(2),foc(3),'k','fill');

% % image plane
add=600;
Ix=uo1*-f/au1+add;
Iy=vo1*-f/av1+add; 
im_cor1=[-Ix,-Iy,f];
im_cor2=[Ix,-Iy,f];
im_cor3=[Ix,Iy,f];
im_cor4=[-Ix,Iy,f];
patch([im_cor1(1) im_cor2(1) im_cor3(1) im_cor4(1)],[im_cor1(2) im_cor2(2) im_cor3(2) im_cor4(2)],[im_cor1(3) im_cor2(3) im_cor3(3) im_cor4(3)],'c');

% % 2D points
points1_2d=V;
points1_2d(1,:)=f*points1_2d(1,:)./points1_2d(3,:);
points1_2d(2,:)=f*points1_2d(2,:)./points1_2d(3,:);
points1_2d(3,:)=f;
scatter3(points1_2d(1,:),points1_2d(2,:),points1_2d(3,:));
for i=1:1   %size(V,2)  , please uncomment "size(V,2)" to draw all the lines that connect the 3D points to their projections on the first image plane
line([0 V(1,i)],[0 V(2,i)],[0 V(3,i)],'color','r');
end

% epipole geometry for the first image plane
L=F_least'*points2;
x1=700;
x2=-700;
m2=zeros(size(V,2),1);
d2=zeros(size(V,2),1);
y1=zeros(size(V,2),1);
y2=zeros(size(V,2),1);
for i=1:size(V,2)
m2(i)=-L(1,i)/L(2,i);
d2(i)=-L(3,i)/L(2,i);
y1(i)=m2(i)*x1+d2(i);
y2(i)=m2(i)*x2+d2(i);
end
epip_x1=(-x1+uo1)*-f/au1;epip_x1=ones(size(V,2),1)*epip_x1;
epip_x2=(-x2+uo1)*-f/au1;epip_x2=ones(size(V,2),1)*epip_x2;
epip_y1=(-y1+vo1)*-f/av1;
epip_y2=(-y2+vo1)*-f/av1;
for i=1:size(V,2)
line([epip_x1(i) epip_x2(i)],[epip_y1(i) epip_y2(i)],[f f],'color','k');
end

% projecting the epipole
Oc=p2*[0;0;0;1];
e(1,1)=f*Oc(1,1)/Oc(3,1);
e(2,1)=f*Oc(2,1)/Oc(3,1);
scatter3(e(1,1),e(2,1),f,'fill','k');


% drawing pi plane
for i=1:1  %size(V,2) , please uncomment "size(V,2)" to draw all PI planes
nor=cross(Oc,V(1:3,i));
nor=nor./norm(nor);
plane_x=2000;plane_y=2000;
planez1=(-plane_x*nor(1)-0*nor(2))/nor(3);
planez2=(plane_x*nor(1)-0*nor(2))/nor(3);
planez3=(0*nor(1)-plane_y*nor(2))/nor(3);
planez4=(0*nor(1)+plane_y*nor(2))/nor(3);
patch([plane_x;0;-plane_x;0],[0;plane_y;0;-plane_y],[planez1;planez3;planez2;planez4],'w');
end


% %  %  %     camera 2 drawing

% computing the origin and standard unit vectors of the second camera's coordinates system
Oc=p2*[0;0;0;1];
Xc=p2*[W;0;0;1];
Yc=p2*[0;W;0;1];
Zc=p2*[0;0;W;1];

% camera x axis
line([Oc(1) Xc(1)],[Oc(2) Xc(2)],[Oc(3) Xc(3)]);hold on;
text(Xc(1),Xc(2),Xc(3),'X axis');

% camera y axis
line([Oc(1) Yc(1)],[Oc(2) Yc(2)],[Oc(3) Yc(3)]);
text(Yc(1),Yc(2),Yc(3),'Y axis');

% camera z axis
line([Oc(1) Zc(1)],[Oc(2) Zc(2)],[Oc(3) Zc(3)]);
text(Zc(1),Zc(2),Zc(3),'Z axis');
text(Oc(1),Oc(2),Oc(3),'Oc2');

% image plane
Ix=uo2*-f/au2+add; % k=f/au
Iy=vo2*-f/av2+add; % k=f/av
im_cor1=[-Ix,-Iy,f];im_cor1=p2*[im_cor1,1]';
im_cor2=[Ix,-Iy,f];im_cor2=p2*[im_cor2,1]';
im_cor3=[Ix,Iy,f];im_cor3=p2*[im_cor3,1]';
im_cor4=[-Ix,Iy,f];im_cor4=p2*[im_cor4,1]';
patch([im_cor1(1) im_cor2(1) im_cor3(1) im_cor4(1)],[im_cor1(2) im_cor2(2) im_cor3(2) im_cor4(2)],[im_cor1(3) im_cor2(3) im_cor3(3) im_cor4(3)],'g');

% % 2d points
points2_2d=ext2*V; % referenced with respect to camera 2
points2_2d(1,:)=f*points2_2d(1,:)./points2_2d(3,:);
points2_2d(2,:)=f*points2_2d(2,:)./points2_2d(3,:);
points2_2d(3,:)=f;
points2_2d=p2*points2_2d;
scatter3(points2_2d(1,:),points2_2d(2,:),points2_2d(3,:),'r');
for i=1:1   %size(V,2)  , please uncomment "size(V,2)" to draw all the lines that connect the 3D points to their projections on the second image plane
line([Oc(1,1) V(1,i)],[Oc(2,1) V(2,i)],[Oc(3,1) V(3,i)],'color','r');
end
line([Oc(1,1) 0],[Oc(2,1) 0],[Oc(3,1) 0],'color','r');
% epipole geometry for the first image plane
L_prime=F_least*points1;
x1=-700;
x2=700;
m=zeros(size(V,2),1);
d=zeros(size(V,2),1);
y1=zeros(size(V,2),1);
y2=zeros(size(V,2),1);
for i=1:size(V,2)
m(i)=-L_prime(1,i)/L_prime(2,i);
d(i)=-L_prime(3,i)/L_prime(2,i);
y1(i)=m(i)*x1+d(i);
y2(i)=m(i)*x2+d(i);
end
epip_x1=(-x1+uo2)*-f/au2;epip_x1=ones(size(V,2),1)*epip_x1;
epip_x2=(-x2+uo2)*-f/au2;epip_x2=ones(size(V,2),1)*epip_x2;
epip_y1=(-y1+vo2)*-f/av2;
epip_y2=(-y2+vo2)*-f/av2;
epip_xy1=p2*[epip_x1';epip_y1';f*ones(1,size(V,2));ones(1,size(V,2))];
epip_xy2=p2*[epip_x2';epip_y2';f*ones(1,size(V,2));ones(1,size(V,2))];
for i=1:size(V,2)
line([epip_xy1(1,i) epip_xy2(1,i)],[epip_xy1(2,i) epip_xy2(2,i)],[epip_xy1(3,i) epip_xy2(3,i)]);
end
% % projecting the epipole
[~,~,v_mat] = svd(F');
e_prime_Fund=v_mat(1:2,end)/v_mat(end,end);
e_prime_Fund(1)=(-e_prime_Fund(1)+uo2)*-f/au2;
e_prime_Fund(2)=(-e_prime_Fund(2)+vo2)*-f/av2;
e_prime_Fund=p2*[e_prime_Fund(1);e_prime_Fund(2);f;1];
scatter3(e_prime_Fund(1),e_prime_Fund(2),e_prime_Fund(3),'fill','k');

title('The entire simulated environment');

