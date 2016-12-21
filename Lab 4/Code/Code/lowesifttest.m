% Gourab Ghosh Roy and Mohammad Rami Koujan

% Script to test the performance of David Lowe's version of SIFT
% Place this file in the folder siftDemoV4

clc
clear all
close all

% Projective Transformation

addpath([cd '/../Sequence1Homographies']);         
load('Sequence1Homographies.mat')
numImages = 16;           % Specify number of images
axisx = 1:1:16;       % Specify values for plotting

figure(); title('Sensitivity to projective transformation');
xlabel('Projective Transformation'); ylabel('Correctly matched (%)'); hold on         % Configure plot

threshold = 2.0;         % threshold to be chosen for matching

noise = ['a','b','c','d'];
for k = 1:4             % for each noise type
    correctMatches = []; 
    for i = 1:numImages    % for each image
        testfilename = strcat('Image_0',num2str(i),noise(k),'.png');    % Get testfilename
        [num,pointMatch] = match('Image_00a.png',testfilename,'false');  % Use a modified version of Lowe's match function
        
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
    switch k
        case 1
            plot(axisx,correctMatches,'-bo'), hold on;
        case 2
            plot(axisx,correctMatches,'-gd'); hold on;
        case 3
            plot(axisx,correctMatches,'-rs'); hold on
        case 4
            plot(axisx,correctMatches,'-k^'); hold on
    end
end

legend('No Noise','With Noise Sigma = 3','With Noise Sigma = 6','With Noise Sigma = 18','Location','southwestoutside'); 
hold off


% Zoom

% Load images and homography accordingly
addpath([cd '/../Sequence2Homographies']);         
load('Sequence2Homographies.mat')
numImages = 9;           % Specify number of images
axisx = 110:5:150;       % Specify values for plotting

figure(); title('Sensitivity to change in camera zoom');
xlabel('Zoom (%)'); ylabel('Correctly matched (%)'); hold on         % Configure plot

threshold = 1.5;         % threshold to be chosen for matching

noise = ['a','b','c','d'];
for k = 1:4             % for each noise type
    correctMatches = [];
    for i = 1:numImages           % for each image
        testfilename = strcat('Image_0',num2str(i),noise(k),'.png');    % Get testfilename
        [num,pointMatch] = match('Image_00a.png',testfilename,'false');  % Use a modified version of Lowe's match function
        
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
    switch k
        case 1
            plot(axisx,correctMatches,'-bo'), hold on;
        case 2
            plot(axisx,correctMatches,'-gd'); hold on;
        case 3
            plot(axisx,correctMatches,'-rs'); hold on
        case 4
            plot(axisx,correctMatches,'-k^'); hold on
    end
end

legend('No Noise','With Noise Sigma = 3','With Noise Sigma = 6','With Noise Sigma = 18','Location','southwestoutside');
hold off


% Rotation

% Load images and homography accordingly
addpath([cd '/../Sequence3Homographies']);         
load('Sequence3Homographies.mat')
numImages = 18;           % Specify number of images
axisx = [-45:5:-5,5:5:45];       % Specify values for plotting

figure(); title('Sensitivity to change in camera rotation');
xlabel('Rotation (degrees)'); ylabel('Correctly matched (%)'); hold on         % Configure plot

threshold = 1.5;         % threshold to be chosen for matching

noise = ['a','b','c','d'];
for k = 1:4             % for each noise type
    correctMatches = [];
    for i = 1:numImages        % for each image
        testfilename = strcat('Image_0',num2str(i),noise(k),'.png');    % Get testfilename
        [num,pointMatch] = match('Image_00a.png',testfilename,'false');  % Use a modified version of Lowe's match function
        
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
    switch k
        case 1
            plot(axisx,correctMatches,'-bo'), hold on;
        case 2
            plot(axisx,correctMatches,'-gd'); hold on;
        case 3
            plot(axisx,correctMatches,'-rs'); hold on
        case 4
            plot(axisx,correctMatches,'-k^'); hold on
    end
end

legend('No Noise','With Noise Sigma = 3','With Noise Sigma = 6','With Noise Sigma = 18','Location','southwestoutside'); 
hold off

