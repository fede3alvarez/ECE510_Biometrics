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

% Separation of local miniman when inspecting Intensity
min_dist = 10;

% Window Parameters:
%   window_size = x_im / window_factor;
%   1st window start = 1;
%   1st window end = 1st window start + window_size;
%   n-th window start = (n-1)th window end - window_overlap;
%   n-th window end = n-th window start + window_size;
window_factor = 10;
window_overlap = 0;

c = 1;

for m = 1:size(Available_images,1)

    %---------------------------------------
    % Step 0: Load and Plot Image
    %---------------------------------------

    % In this case, it is in Grayscale
    current_image = Available_images(m,:);
    [Image_A, Map_A] = imread(current_image);

    [y_im, x_im] = size(Image_A);
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

    x_gauss = [];
    y_gauss = [];

    % The "user-assisted initialization" referenced in the paper
    %   will be interpreted to be the lower intensity detected
    x_int = [];
    y_int = [];

    % From paper:
    %    pvi_Di_yi = pvi_Di * pvi_yi
    %    pvi_Di = c * {1 - Di / (max_k * Dk)}
    %    pvi_yi = {1 / [sigma*sqrt(2*pi)]} * exp{-(yi-y_int)^2/(sigma^2)}
    pvi_yi = [];
    pvi_Di = [];
    pvi_Di_yi = [];

    % Iterate through "windows" until the image is swept
    while (x_im > window_end)

        % For every window sweep,
        %   Iterate over the y axxis, and get the avrg intensity
        %for y_curr = round(y_im/4):round(y_im*3/4)
        for y_curr = round(y_im/6):round(y_im*5/6)
            y_int = [y_int, y_curr];
            x_int = [x_int, mean(...
                             Image_A(y_curr,[window_start:window_end])...
                            )];
        end
    
        % Smooth Intensity
        x_int = smoothn(x_int);

        % Plot Image
        f_figure = f_figure + 1;
        figure(f_figure)
        plot(x_int, y_int)

        [min_val, min_idx] = min(x_int);
        [max_val, max_idx] = max(x_int);
        
        x = [x, mean(window_start,window_end)];
        y = [y, y_int(min_idx)];

    
        %--------------------------------------------
        % Step 2: Calculate Pvi(yi)
        %   Use of Exponentials
        %--------------------------------------------
        
        % Finding sigma
        % https://ned.ipac.caltech.edu/level5/Leo/Stats2_3.html
        % We will assume that the distance from the max to 
        %   the beginning or end (which ever is shorter) covers 
        %   99.7% of the Gaussian and is equal to 6 sigma
        sigma = min(max_idx,(window_size - max_idx))/6;
        
        % Finding all local minima
        int_min = islocalmin(x_int,'MinSeparation',min_dist);
        x_int_min = x_int(int_min);
        y_int_min = y_int(int_min);

        % We'll assume that y_hat is always the minimun closer to the
        % middle
        y_hat = y_im / 2;
        maxK = sum(cumtrapz(x_int,y_int));

        pvi_yi = [];
        pvi_Di = [];
        pvi_Di_yi = [];

        % Iterate through each local minima
        for n = 1 : length(x_int_min)
            x_i = x_int_min(n);
            y_i = y_int_min(n);
    
            % Finding pvi_yi
            % From paper:
            % pvi_yi = {1 / [sigma*sqrt(2*pi)])} * 
            %               exp{-(yi-y_hat)^2/(sigma^2)}
            pvi_yi_curr = (1/(sigma*sqrt(2*pi)))* ...
                            exp(-(y_i-y_hat)^2 / (sigma^2));
            pvi_yi = [pvi_yi , pvi_yi_curr ];
            
            % Finding pvi_Di
            % From paper:
            %    pvi_Di = c * {1 - Di / (max_k * Dk)}
            D_i_curr = (window_end - window_start) - x_i;
            pvi_Di_curr = c * (1 - D_i_curr) * maxK;
            pvi_Di = [pvi_Di, pvi_Di_curr];

            % From paper:
            %    pvi_Di_yi = pvi_Di * pvi_yi
            pvi_Di_yi_curr = pvi_yi_curr * pvi_Di_curr;
            pvi_Di_yi = [pvi_Di_yi, pvi_Di_yi_curr];

        end

        [gauss_val, gauss_idx] = max(pvi_Di_yi);
    
        x_gauss = [x_gauss, x_int_min(gauss_idx)];
        y_gauss = [y_gauss, y_int_min(gauss_idx)];

        %--------------------------------------------
        % Update parameters for next cycle
        %--------------------------------------------

        window_start = window_end - window_overlap;
        window_end = window_start + window_size;
        x_int = [];
        y_int = [];

        if (window_end >= x_im)
            window_end = x_im;
        end

    %--------------------------------------------
    % Step 3: Sweep the Window
    %--------------------------------------------
    
    end

    % Plot Image
    figure(1)
    hold on
    plot(x,y,'g*')
    hold on
    plot(x_gauss,y_gauss,'m*')


    %--------------------------------------------
    % Step 4: 2nd Degree Polynomial
    %--------------------------------------------
    
    coeff_valley = polyfit(x,y,2);
    x_valley = 1:x_im;
    y_valley = polyval(coeff_valley, x_valley);

    % Plot Image
    figure(1)
    hold on
    plot(x_valley, y_valley, 'b-')

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
