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
 Available_images = ['iris1.png'
                    'iris2.png'
                    'iris3.png'
                    'iris4.png'
                    'iris5.png'
                    'iris6.png'
                    ];


f_figure = 1;

for m = 1:size(Available_images,1)
    current_image = Available_images(m,:);
    [Image_A, Map_A] = imread(current_image);

    size(Image_A)
    size(Map_A)

    % Convert to Grayscale
    % Gray_Map_A = rgb2gray(Image_A);
    % size(Gray_Map_A)
    % Filter Image
    %Gray_Map_A = imgaussfilt(Gray_Map_A,5);
    
    % Plot Image
    % figure(f_figure)

end

% Fix Matlab axis so getpts is accurate
    % hax = axes('Parent', figure(f_figure));
    % axis(hax,'manual');
    % imshow(Image_A);
    % title(current_image);
        
    % % Get user to select Point
    %[x,y] = getpts;

    % % Save data
    % save_file_selections = strcat(current_image(1:5),'_data')
    % save(save_file_selections,'x','y')

    %% Load saved data
    % load_selections = strcat(current_image(1:11),'_data.mat');
    % selections = importdata(load_selections);
    % x = selections.x;
    % y = selections.y;
    % hold on

% 
%     f_figure = f_figure + 1;
% end         % All images Iteration
