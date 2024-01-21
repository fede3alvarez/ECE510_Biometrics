close all
clear all;
clc;

% Code will be available at
% https://github.com/fede3alvarez/ECE510_Biometrics

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Read Image (*.jpg) and massage it       %
%             to make matlab happy...            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%/home/fico/Repos/ECE510_Biometrics/HW_2/Images/00_Hand.jpg\

% Note: jpeg images are read directly as RGB  matrices
%    https://www.mathworks.com/matlabcentral/answers/23119-jpg-to-rgb-image
[Image_A, Map_A] = imread('00_Hand.jpg');
Image = double(Image_A);

% Separete data in colors
[rows, cols, colors]  = size(Image);
red = Image_A(:, :, 1);
green = Image_A(:, :, 2);
blue = Image_A(:, :, 3);

% Rearrange data for kmeans clustering
data = double([red(:), green(:), blue(:)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              K-Means Clustering                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In this section we setup the numbers of K cluster to be used

% TO-DO: Implement this as a for-loop (2,3,5)
k = 2;

% Setting random k initial conditions
init_cond = []

for i = 1:k 
    r = randi(size(data,1),1);
    init_cond = [init_cond; data(r,:)];
end

[Hand, Background] = kmeans(data, k,'start', init_cond)

%indexes = kmeans(data, 2);
%class1 = reshape(indexes == 1, rows, cols);



% Init_Values = [r0, g0, b0;r1, g1, b1]
% Init_Values = [r0, g0, b0;r1, g1, b1;r2, g2, b2]
% [Hand, Background] = kmeans(Image, 3,'start', Init_Values)
% % WORKING:
% %[Hand, Background] = kmeans(Image(:,:,1), 2)

%class1 = reshape(indexes == 1, rows, cols);
% 
% % Plots for debugging
% subplot(2, 2, 1);
% imshow(Image_A);
% title('Original Color Image');
% subplot(2, 2, 2);
% imshow(red);
% title('Red Image');
% subplot(2, 2, 3);
% imshow(green);
% title('Green Image');
% subplot(2, 2, 4);
% imshow(blue);
% title('Blue Image');
