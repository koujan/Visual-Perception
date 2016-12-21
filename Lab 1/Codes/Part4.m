% part 4
clear all;
clc
window_size=11; % suppression window size 
dirc=dir(strcat(pwd,'\*.png')); % the images should be in the same working directory 
for i=1:length(dirc)
    im=imread(dirc(i).name);
    [E]=compute_matrix_E(im);
    figure;
    imshow(mat2gray(E));hold on;
    title(horzcat('non-maximal suppression for Matrix E of ',dirc(i).name,' image'));
    features(81)=struct('p_x',[],'p_y',[]);
    E_modified=padarray(E,[floor(window_size/2) floor(window_size/2)],min(min(E))); % adding borders with the minimum value of E to the original E image to make sure that opening window around the maximum pixel will not go outside the image if the same pixel is near the borders
    for i=1:81
        [val,ind]=max(E_modified);
        [val2,ind2]=max(val);
        features(i).p_x=ind(ind2)-floor(window_size/2);
        features(i).p_y=ind2-floor(window_size/2);
        E_modified(ind(ind2),ind2)=min(min(E_modified));
        plot(features(i).p_y, features(i).p_x, 'r+');
        text(features(i).p_y, features(i).p_x,horzcat('   ',num2str(i)),'color','red');
        E_modified(ind(ind2)-floor(window_size/2):ind(ind2)+floor(window_size/2),ind2-floor(window_size/2):ind2+floor(window_size/2))=min(min(E));
    end
    hold off;
end

% the same procedure but for matrix R
clear all;clc
window_size=25; % suppression window size
dirc=dir(strcat(pwd,'\*.png')); % the images should be in the same working directory
for i=1:length(dirc)
    im=imread(dirc(i).name);
    [R]=compute_matrix_R(im);
    figure;
    imshow(mat2gray(R));hold on;
    title(horzcat('non-maximal suppression for Matrix R of ',dirc(i).name,' image'));
    features(81)=struct('p_x',[],'p_y',[]);
    R_modified=padarray(R,[floor(window_size/2) floor(window_size/2)],min(min(R))); % adding borders with the minimum value of E to the original E image to make sure that opening window around the maximum pixel will not go outside the image if the same pixel is near the borders
    for i=1:81
        [val,ind]=max(R_modified);
        [val2,ind2]=max(val);
        features(i).p_x=ind(ind2)-floor(window_size/2);
        features(i).p_y=ind2-floor(window_size/2);
        R_modified(ind(ind2),ind2)=min(min(R_modified));
        plot(features(i).p_y, features(i).p_x, 'r+');
        text(features(i).p_y, features(i).p_x,horzcat('   ',num2str(i)),'color','red');
        R_modified(ind(ind2)-floor(window_size/2):ind(ind2)+floor(window_size/2),ind2-floor(window_size/2):ind2+floor(window_size/2))=min(min(R));
    end
    hold off;
end