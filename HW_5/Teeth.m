close all
clear all;
clc;

% Code will be available at
% https://github.com/fede3alvarez/ECE510_Biometrics
% Homework 5 - Teeth


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Images
 Available_images = ['TeethSample.png'];

% All images Iteration
% This Homework has a single image...
f_figure = 1;

% Window Parameters:
%   window_size = x_im / window_factor;
%   1st window start = 1;
%   1st window end = 1st window start + window_size;
%   n-th window start = (n-1)th window end - window_overlap;
%   n-th window end = n-th window start + window_size;
window_factor = 15;
window_overlap = 0;

sigma = 1;

for m = 1:size(Available_images,1)

    %---------------------------------------
    % Step 0: Load and Plot Image
    %---------------------------------------

    % In this case, it is in Grayscale
    current_image = Available_images(m,:);
    [Image_A, Map_A] = imread(current_image);

    [y_hat_im, x_im] = size(Image_A);
    [y_map, x_map] = size(Map_A);
  
    % Plot Image
    figure(f_figure)
    %Fix Matlab axis if getpts is needed, it is accurate
    % hax = axes('Parent', figure(f_figure));
    % axis(hax,'manual');
    imshow(Image_A);
    title(current_image);

    % f_figure = f_figure + 1;


    %--------------------------------------------
    % Step 1: Define Sliding Window
    %--------------------------------------------

    % Define Window
    window_size = round(x_im / window_factor);
    window_start = 1;
    window_end = window_start + window_size - window_overlap;

    % Initialize arrays to collect windows analysis results
    x = [];
    y = [];

    % The "user-assisted initialization" referenced in the paper
    %   will be interpreted to be the lower intensity detected
    int_i = [];
    y_hat_i = [];

    % From paper:
    %    pvi_Di_yi = pvi_Di * pvi_yi
    %    pvi_Di = c * {1 - Di / (max_k * Dk)}
    %    pvi_yi = {1 / sqrt(2*pi*sigma)} * exp{-(yi-y_hat_i)^2/(sigma^2)}
    pvi_yi = [];
    pvi_Di = [];
    pvi_Di_yi = [];

    % Iterate through "windows" until the image is swept
    while (x_im > window_end)

        % For every window sweep,
        %   Iterate over the y axxis, and get the avrg intensity
        for y_curr = round(y_hat_im/4):round(y_hat_im*3/4)
            y_hat_i = [y_hat_i, y_curr];
            int_i = [int_i, mean(...
                             Image_A(y_curr,[window_start:window_end])...
                            )];
        end
    
        % Plot Image
        f_figure = f_figure + 1;
        figure(f_figure)
        plot(int_i, y_hat_i)

        [min_val, min_idx] = min(int_i);
        x = [x, mean(window_start,window_end)];
        y = [y, y_hat_i(min_idx)];

    
        %--------------------------------------------
        % Step 2: Calculate Pvi(yi)
        %   Use of Exponentials
        %--------------------------------------------
        
        % From paper:
        % pvi_yi = {1 / sqrt(2*pi*sigma)} * exp{-(yi-y_hat_i)^2/(sigma^2)}
        
        pvi_yi_curr = 0;
        pvi_yi = [pvi_yi , pvi_yi_curr ];
            

        %--------------------------------------------
        % Update parameters for next cycle
        %--------------------------------------------

        window_start = window_end - window_overlap;
        window_end = window_start + window_size;
        int_i = [];
        y_hat_i = [];

        if (window_end >= x_im)
            window_end = x_im;
        end

    end

    % Plot Image
    figure(1)
    hold on
    plot(x,y,'g*')


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
