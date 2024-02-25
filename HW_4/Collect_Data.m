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
 Available_images = [%'iris1.jpg'
%                     'iris2.jpg'
%                     'iris3.jpg'
%                     'iris4.jpg'
%                     'iris5.jpg'
                    'iris6.jpg'
                    ];


f_figure = 1;

for m = 1:size(Available_images,1)
    current_image = Available_images(m,:);
    [Image_A, Map_A] = imread(current_image);

    % Convert to Grayscale
    %Gray_Map_A = rgb2gray(Image_A);
    
    % Filter Image
    %Gray_Map_A = imgaussfilt(Gray_Map_A,5);
    
    % Plot Image
    figure(f_figure)

    % Fix Matlab axis so getpts is accurate
    hax = axes('Parent', figure(f_figure));
    axis(hax,'manual');
    imshow(Image_A);
    title(current_image);
        
    % % Get user to select Point
    [x,y] = getpts;

    % % Save data
    save_file_selections = strcat(current_image(1:5),'_data')
    save(save_file_selections,'x','y')

    %% Load saved data
    % load_selections = strcat(current_image(1:11),'_data.mat');
    % selections = importdata(load_selections);
    % x = selections.x;
    % y = selections.y;
    % hold on


    f_figure = f_figure + 1;
end         % All images Iteration
