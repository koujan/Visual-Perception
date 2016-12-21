% part 5
clear all;clc;
window_size=11;
dirc=dir(strcat(pwd,'\*.png')); % the images should be in the same working directory
for j=1:length(dirc)
    im=imread(dirc(j).name);
    [E]=compute_matrix_E(im);
    figure;imshow(im);hold on;
    title(horzcat('Subpixel accuracy method for ',dirc(j).name,' image based on E matrix'));
    features(81)=struct('p_x',[],'p_y',[]);
    E_modified=padarray(E,[floor(window_size/2) floor(window_size/2)],min(min(E))); % adding borders with the minimum value of E to the original E image to make sure that opening window around the maximum pixel will not go outside the image if the same pixel is near the borders
    for i=1:81
        [val,ind]=max(E_modified);
        [val2,ind2]=max(val);
        features(i).p_x=ind(ind2)-floor(window_size/2); % removing the effect of the padding
        features(i).p_y=ind2-floor(window_size/2);
        %E_modified(ind(ind2),ind2)=min(min(E_modified));% in order not to be selected again
        E_modified(ind(ind2)-floor(window_size/2):ind(ind2)+floor(window_size/2),ind2-floor(window_size/2):ind2+floor(window_size/2))=min(min(E)); % avoiding selection of neighbours again
       
        % actuall start of part 5
  
        f=@(x,y)[x.^2; y.^2; x; y; x.*y; x-x+1];u=0;
        A=f([-1 -1 -1 0 0 0 1 1 1],repmat([-1:1],[1 3])); A=A'; % building the A matrix
        E_modified2=padarray(E,[floor(window_size/2) floor(window_size/2)]);
        b=E_modified2(features(i).p_x+floor(window_size/2)-1:features(i).p_x+floor(window_size/2)+1,features(i).p_y+floor(window_size/2)-1:features(i).p_y+floor(window_size/2)+1);
        b=b(:);
        X=A\b; % solving the system
        % calculating the maximum of the parabola by brute force approach
        Fxy=@(x,y,a,b,c,d,e,f) a*x.^2+b*y.^2+c*x+d*y+e.*x.*y+f;
        [x,y] = meshgrid(-1:.01:1,-1:.01:1);
        [val,ind]=max(Fxy(x,y,X(1),X(2),X(3),X(4),X(5),X(6)));
        [val2,ind2]=max(val);
        x_1=ind(ind2);
        x_2=ind2;
        plot(features(i).p_y+x(x_1,x_2),features(i).p_x+y(x_1,x_2) ,'r+');
    end
end


% the same procedure but for matrix R

clear all;clc;
window_size=11;
dirc=dir(strcat(pwd,'\*.png')); % the images should be in the same working directory
for j=1:length(dirc)
    im=imread(dirc(j).name);
    [R]=compute_matrix_R(im);
    figure;imshow(im);hold on;
    title(horzcat('Subpixel accuracy method for ',dirc(j).name,' image based on R matrix'));
    features(81)=struct('p_x',[],'p_y',[]);
    R_modified=padarray(R,[floor(window_size/2) floor(window_size/2)],min(min(R))); % adding borders with the minimum value of E to the original E image to make sure that opening window around the maximum pixel will not go outside the image if the same pixel is near the borders
    for i=1:81
        [val,ind]=max(R_modified);
        [val2,ind2]=max(val);
        features(i).p_x=ind(ind2)-floor(window_size/2); % removing the effect of the padding
        features(i).p_y=ind2-floor(window_size/2);
        %E_modified(ind(ind2),ind2)=min(min(E_modified));% in order not to be selected again
        R_modified(ind(ind2)-floor(window_size/2):ind(ind2)+floor(window_size/2),ind2-floor(window_size/2):ind2+floor(window_size/2))=min(min(R)); % avoiding selection of neighbours again
       
        % actuall start of part 5
  
        f=@(x,y)[x.^2; y.^2; x; y; x.*y; x-x+1];u=0;
        A=f([-1 -1 -1 0 0 0 1 1 1],repmat([-1:1],[1 3])); A=A'; % building the A matrix
        R_modified2=padarray(R,[floor(window_size/2) floor(window_size/2)]);
        b=R_modified2(features(i).p_x+floor(window_size/2)-1:features(i).p_x+floor(window_size/2)+1,features(i).p_y+floor(window_size/2)-1:features(i).p_y+floor(window_size/2)+1);
        b=b(:);
        X=A\b; % solving the system
        % calculating the maximum of the parabola by brute force approach
        Fxy=@(x,y,a,b,c,d,e,f) a*x.^2+b*y.^2+c*x+d*y+e.*x.*y+f;
        [x,y] = meshgrid(-1:.01:1,-1:.01:1);
        [val,ind]=max(Fxy(x,y,X(1),X(2),X(3),X(4),X(5),X(6)));
        [val2,ind2]=max(val);
        x_1=ind(ind2);
        x_2=ind2;
        plot(features(i).p_y+x(x_1,x_2),features(i).p_x+y(x_1,x_2) ,'r+');
    end
end


