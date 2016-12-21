function [E]=compute_matrix_E(im)

%im = imread('chessboard04.png');
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

%%% part 1
m1=conv2(Ix2,ones(3,3),'same');
m2=conv2(Iy2,ones(3,3),'same');
m3=conv2(Ixy,ones(3,3),'same');
E=zeros(size(im));
for i=1:size(im,1)
    for j=1:size(im,2)
        M=[m1(i,j) m3(i,j); m3(i,j) m2(i,j)];
        [V,D] = eig(M);
        E(i,j)=min(D(1,1),D(2,2));
    end
end

end