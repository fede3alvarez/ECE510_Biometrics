close all
clear all;
clc;

% Code will be available at
% https://github.com/fede3alvarez/ECE510_Biometrics
% Homework 4 - Iris

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Data Collection                 %
%            (Ask Users to Select               %
%            Eye / Iris Center and              %
%             saves it to a file)               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Images
 Available_images = ['TeethSample.png'];

% All images Iteration
% This Homework has a single image...
f_figure = 1;

% Window Parameters:
%   window_size = col_im / window_factor;
%   1st window start = 0;
%   1st window end = 1st window start + window_size;
%   n-th window start = (n-1)th window end - window_overlap;
%   n-th window end = n-th window start + window_size;
window_factor = 10;
window_overlap = 0;

sigma = 1;

for m = 1:size(Available_images,1)

    %---------------------------------------
    % Step 0: Load and Plot Image
    %---------------------------------------

    % In this case, it is in Grayscale
    current_image = Available_images(m,:);
    [Image_A, Map_A] = imread(current_image);

    [row_im, col_im] = size(Image_A);
    [row_map, col_map] = size(Map_A);
  
    % Plot Image
    figure(f_figure)
    %Fix Matlab axis if getpts is needed, it is accurate
    hax = axes('Parent', figure(f_figure));
    axis(hax,'manual');
    imshow(Image_A);
    title(current_image);

    %f_figure = f_figure + 1;


    %--------------------------------------------
    % Step 1: Define Sliding Window
    %--------------------------------------------

    %--------------------------------------------
    % Step 2: Calculate Pvi(yi)
    %--------------------------------------------

    %--------------------------------------------
    % Step 3: Sweep the Window
    %--------------------------------------------
    
    %--------------------------------------------
    % Step 4: 2nd Degree Polynomial
    %--------------------------------------------
    
    %--------------------------------------------
    % Step 5: Repeat Steps 2-4 Horizontally
    %--------------------------------------------

    %--------------------------------------------
    % Step 6: Repeat Steps for Upper & Lower
    %              Ranges Separetely
    %--------------------------------------------

    %--------------------------------------------
    % Step 7: Teeth Separation
    %--------------------------------------------

end         % All images Iteration
