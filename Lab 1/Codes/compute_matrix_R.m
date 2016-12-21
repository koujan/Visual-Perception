function [R]=compute_matrix_R(im)

im_orig = im;
if size(im,3)>1 im=rgb2gray(im); end
% Derivative masks
dx = [-1 0 1;   
-1 0 1;
-1 0 1];
dy = dx';
% Image derivatives
Ix = conv2(double(im), dx, 'same');
Iy = conv2(double(im), dy, 'same');
sigma=2;
% Generate Gaussian filter of size 9x9 and std. dev. sigma.
g = fspecial('gaussian',9, sigma);
% Smoothed squared image derivatives
Ix2 = conv2(Ix.^2, g, 'same');
Iy2 = conv2(Iy.^2, g, 'same');
Ixy = conv2(Ix.*Iy, g, 'same');

%%% part 2
m1=conv2(Ix2,ones(3,3),'same');
m2=conv2(Iy2,ones(3,3),'same');
m3=conv2(Ixy,ones(3,3),'same');
R=zeros(size(im));
for i=1:size(im,1)
    for j=1:size(im,2)
        M=[m1(i,j) m3(i,j); m3(i,j) m2(i,j)];
        R(i,j)=(M(1,1)*M(2,2)-M(1,2)*M(2,1))-0.04*(M(1,1)+M(2,2))^2;
        %R(i,j)=det(M)-0.04*trace(M);
    end
end

end