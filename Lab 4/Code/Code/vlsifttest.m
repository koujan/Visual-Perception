% SIFT Test based on VLFeat

clc
clear all
close all
% Run VL Setup, will change accordingly on a different computer
run C:\Vlfeat\vlfeat-0.9.20-bin.tar\vlfeat-0.9.20-bin\vlfeat-0.9.20\toolbox\vl_setup


% Projective Transformation

% Load images and homography accordingly
addpath([cd '/Sequence1Homographies']);  % ensuring directories are at this level
load('Sequence1Homographies.mat')
numImages = 16;         % Specify number of images
axisx = 1:1:16;     % Specify values for plotting

threshold = 2.5;        % threshold to be chosen for matching

magnif = [3.0,2.2];     % Values of magnification factor to change descriptor size

noise = ['a','b','c','d'];    % different noise levels
disp('Projective Transformation')
for k = 1:4       % for each noise level
    figure(); title('Sensitivity to projective transformation');
    xlabel('Projective Transformation'); ylabel('Correctly matched (%)'); hold on
    for l = 1:2    % for each descriptor size
        correctMatches = [];
        for i = 1:numImages  % for each image 
            testfilename = strcat('Image_0',num2str(i),noise(k),'.png');    % Get testfilename
            [num,pointMatch] = vlmatch('Image_00a.png',testfilename,magnif(l));  % Use a modified version of Lowe's match function
            
            correctMatch = 0;
            for j = 1:num
                p_00 = [pointMatch(j,1) pointMatch(j,2) 1];
                shift=200;
                % Get corresponding point from homography
                if (i<9)
                    p_01 = Sequence1Homographies(i).H *[1,0,shift;0,1,-1;0,0,1] * p_00';
                    p_01 = p_01/(p_01(3));
                    p_01=[1,0,-shift;0,1,1;0,0,1]*p_01;
                else
                    p_01 = Sequence1Homographies(i).H*[1,0,0;0,1,shift;0,0,1]*p_00';
                    p_01 = p_01/(p_01(3));
                    p_01=[1,0,0;0,1,-shift;0,0,1]*p_01;
                end
                if (abs(p_01(1)-pointMatch(j,3)) <= threshold) && (abs(p_01(2)-pointMatch(j,4)) <= threshold)  %  consider match
                    correctMatch = correctMatch + 1;
                end
            end
            correctMatch = 100*correctMatch/num;        % Get percentage of correct match
            correctMatches = [correctMatches, correctMatch];   % Store in an array for that image
        end
        s = sprintf('Image type = %s',noise(k));
        disp(s)
        s = sprintf('Magnif value = %f',magnif(l));
        disp(s)
        disp('% of Correct matches')
        disp(correctMatches)
        switch l
            case 1
                plot(axisx,correctMatches,'-bo'), hold on;
            case 2
                plot(axisx,correctMatches,'-rs'); hold on;
        end
    end
    legend('16 by 16 descriptor','12 by 12 descriptor','Location','southwest');
    hold off
end



% Zoom

% Load images and homography accordingly
addpath([cd '/Sequence2Homographies']);  % ensuring directories are at this level
load('Sequence2Homographies.mat')
numImages = 9;         % Specify number of images
axisx = 110:5:150;     % Specify values for plotting

threshold = 2.5;        % threshold to be chosen for matching

magnif = [3.0,2.2];     % Values of magnification factor to change descriptor size

noise = ['a','b','c','d'];    % different noise levels
disp('Zoom')
for k = 1:4        % for each noise level
    figure(); title('Sensitivity to change in camera zoom');
    xlabel('Zoom (%)'); ylabel('Correctly matched (%)'); hold on
    for l = 1:2     % for each descriptor size
        correctMatches = [];
        for i = 1:numImages     % for each image 
            testfilename = strcat('Image_0',num2str(i),noise(k),'.png');    % Get testfilename
            [num,pointMatch] = vlmatch('Image_00a.png',testfilename,magnif(l));  % Use a modified version of Lowe's match function
            
            correctMatch = 0;
            for j = 1:num
                p_00 = [pointMatch(j,1) pointMatch(j,2) 1];
                p_01 = Sequence2Homographies(i).H * p_00';          % Get corresponding point from homography
                if (abs(p_01(1)-pointMatch(j,3)) <= threshold) && (abs(p_01(2)-pointMatch(j,4)) <= threshold)  %  consider match
                    correctMatch = correctMatch + 1;
                end
            end
            correctMatch = 100*correctMatch/num;        % Get percentage of correct match
            correctMatches = [correctMatches, correctMatch];   % Store in an array for that image
        end
        s = sprintf('Image type = %s',noise(k));
        disp(s)
        s = sprintf('Magnif value = %f',magnif(l));
        disp(s)
        disp('% of Correct matches')
        disp(correctMatches)
        switch l
            case 1
                plot(axisx,correctMatches,'-bo'), hold on;
            case 2
                plot(axisx,correctMatches,'-rs'); hold on;
        end
    end
    legend('16 by 16 descriptor','12 by 12 descriptor','Location','southwest');
    hold off
end



% Rotation

% Load images and homography accordingly
addpath([cd '/Sequence3Homographies']);  % ensuring directories are at this level
load('Sequence3Homographies.mat')
numImages = 18;          % Specify number of images
axisx = [-45:5:-5,5:5:45];   % Specify values for plotting

threshold = 2.5;              % threshold to be chosen for matching

magnif = [3.0,2.2];           % Values of magnification factor to change descriptor size

noise = ['a','b','c','d'];   % different noise levels
disp('Rotation')
for k = 1:4            % for each noise level
    figure(); title('Sensitivity to change in camera rotation');
    xlabel('Rotation (degrees)'); ylabel('Correctly matched (%)'); hold on
    for l = 1:2        % for each descriptor size
        correctMatches = [];
        for i = 1:numImages         % for each image 
            testfilename = strcat('Image_0',num2str(i),noise(k),'.png');    % Get testfilename
            [num,pointMatch] = vlmatch('Image_00a.png',testfilename,magnif(l));  % Use a modified version of Lowe's match function
            
            correctMatch = 0;
            for j = 1:num
                p_00 = [pointMatch(j,1) pointMatch(j,2) 1];
                p_01 = Sequence3Homographies(i).H * p_00';          % Get corresponding point from homography
                if (abs(p_01(1)-pointMatch(j,3)) <= threshold) && (abs(p_01(2)-pointMatch(j,4)) <= threshold)  %  consider match
                    correctMatch = correctMatch + 1;
                end
            end
            correctMatch = 100*correctMatch/num;        % Get percentage of correct match
            correctMatches = [correctMatches, correctMatch];   % Store in an array for that image
        end
        s = sprintf('Image type = %s',noise(k));
        disp(s)
        s = sprintf('Magnif value = %f',magnif(l));
        disp(s)
        disp('% of Correct matches')
        disp(correctMatches)
        switch l
            case 1
                plot(axisx,correctMatches,'-bo'), hold on;
            case 2
                plot(axisx,correctMatches,'-rs'); hold on;
        end
    end
    legend('16 by 16 descriptor','12 by 12 descriptor','Location','southwest');
    hold off
end



