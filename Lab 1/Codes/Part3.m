% part 3
clear all;
clc
dirc=dir(strcat(pwd,'\*.png')); % the images should be in the same working directory 
for i=1:length(dirc)
    im=imread(dirc(i).name);
    [E]=compute_matrix_E(im);
    subplot(2,2,i);
    imshow(mat2gray(E));hold on;
    features(81)=struct('p_x',[],'p_y',[]); 
    title(horzcat('The 81 most salient points of E matrix of ',dirc(i).name,' image'));
    for c=1:81
        [val,ind]=max(E);
        [val2,ind2]=max(val);
        features(c).p_x=ind(ind2);
        features(c).p_y=ind2;
        E(ind(ind2),ind2)=min(min(E));
        plot(features(c).p_y, features(c).p_x, 'r+');
    end
    hold off;
end

% the same procedure but for matrix R
clear all;clc
figure;
dirc=dir(strcat(pwd,'\*.png')); % the images should be in the same working directory
for i=1:length(dirc) % first two elements of the dirc variable are special characters
    im=imread(dirc(i).name);
    [R]=compute_matrix_R(im);
    subplot(2,2,i);
    imshow(mat2gray(R));hold on;
    title(horzcat('The 81 most salient points of R matrix of  ',dirc(i).name,' image'));
    features(81)=struct('p_x',[],'p_y',[]); 
    for c=1:81
        [val,ind]=max(R);
        [val2,ind2]=max(val);
        features(c).p_x=ind(ind2);
        features(c).p_y=ind2;
        R(ind(ind2),ind2)=min(min(R));
        plot(features(c).p_y, features(c).p_x, 'r+');
    end
    hold off;
   
end