% part 2
clear all;clc
dirc=dir(strcat(pwd,'\*.png')); % the images should be in the same working directory 
for i=1:length(dirc) % first two elements of the dirc variable are special characters
    im=imread(dirc(i).name);
    [R]=compute_matrix_R(im);
    subplot(2,4,i);
    imshow(im);
    title(dirc(i).name);
    subplot(2,4,i+4);
    imshow(mat2gray(R));
    title(horzcat('Matrix R of ',dirc(i).name,' image'));
   
end

