% Part 1
clear all;
clc
dirc=dir(strcat(pwd,'\*.png')); % the images should be in the same working directory

for i=1:length(dirc) % first two elements of the dirc variable are special characters
    im=imread(dirc(i).name);
    [E]=compute_matrix_E(im);
    subplot(2,4,i);
    imshow(im);
    title(dirc(i).name);
    subplot(2,4,i+4);
    imshow(mat2gray(E));
    title(strcat('Matrix E of:',dirc(i).name));
end