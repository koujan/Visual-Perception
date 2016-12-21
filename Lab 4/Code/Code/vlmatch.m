function [num,pointMatch] = vlmatch(image1, image2, magnif)
% function to match SIFT features of 2 input images
% magnif determines descriptor size

% Find SIFT keypoints and match them using VLFeat functions
[f1, des1] = vl_sift(single(imread(image1)),'Magnif', magnif);
[f2, des2] = vl_sift(single(imread(image2)),'Magnif', magnif);
matches = vl_ubcmatch(des1,des2,1.4);

num = size(matches,2);
pointMatch = [];
for i = 1:num
    pointMatch=[pointMatch; [f1(1,matches(1,i)),f1(2,matches(1,i)),f2(1,matches(2,i)),f2(2,matches(2,i))]];  % Save x y values of matches
end

end