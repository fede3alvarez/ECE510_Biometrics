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

%---------------------------------------
% Parameters for Valley calculation
%---------------------------------------

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

% Other Useful parameters
c = 1;
horiz_sample_pts = 150;
numb_of_teeth = 6;


for m = 1:size(Available_images,1)

    %---------------------------------------
    % Step 0: Load and Plot Image
    %---------------------------------------

    % In this case, it is in Grayscale
    current_image = Available_images(m,:);
    [Image_A, Map_A] = imread(current_image);

    Image_A = imgaussfilt(Image_A,3);      
    [y_im, x_im] = size(Image_A);
  
    % Plot Image
    figure(f_figure)
    imshow(Image_A);
    title(current_image);


    %--------------------------------------------
    % Step 1: Define Sliding Window
    %--------------------------------------------

    % Define Window
    window_size = round(x_im / window_factor);
    window_start = round(window_size / 2);
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
        x_curr = mean(window_start,window_end);
        
        % % Plot Image
        % f_figure = f_figure + 1;
        % figure(f_figure)
        % plot(x_int, y_int)

        [min_val, min_idx] = min(x_int);
        [max_val, max_idx] = max(x_int);
        
        x = [x, x_curr];
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

        % We'll assume that y_hat is intensity minimun
        y_hat = y_int(min_idx);
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

        x_gauss = [x_gauss, x_curr];
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

    % % Plot Image
    % figure(1)
    % hold on
    % plot(x_gauss,y_gauss,'m*')

    %--------------------------------------------
    % Step 4: 2nd Degree Polynomial
    %--------------------------------------------
    
    coeff_gauss_valley = polyfit(x_gauss,y_gauss,2);
    x_gauss_valley = 1:x_im;
    y_gauss_valley = polyval(coeff_gauss_valley, x_gauss_valley);

    % Plot Image
    figure(1)
    hold on
    plot(x_gauss_valley, y_gauss_valley, 'm-','LineWidth',2)

    %--------------------------------------------
    % Step 5: Repeat Steps 2-4 Horizontally
    % Step 6: Repeat Steps for Upper & Lower
    %              Ranges Separetely
    %--------------------------------------------

    %--------------------------------------------
    % Step 5.1: Define Sliding Window
    %           for Upper Teeth
    %--------------------------------------------

    % Initialize arrays to collect windows analysis results
    upper_x_curr = [];
    upper_y_curr = [];

    upper_x_teeth_sep = [];
    upper_y_teeth_sep = [];

    % The "user-assisted initialization" referenced in the paper
    %   will be interpreted to be the lower intensity detected
    upper_x_int = [];
    upper_y_int = [];

    % From paper:
    %    pvi_Di_yi = pvi_Di * pvi_yi
    %    pvi_Di = c * {1 - Di / (max_k * Dk)}
    %    pvi_yi = {1 / [sigma*sqrt(2*pi)]} * exp{-(yi-y_int)^2/(sigma^2)}
    upper_pvi_yi = [];
    upper_pvi_Di = [];
    upper_pvi_Di_yi = [];

    %--------------------------------------------
    % Step 5.1: Define Sliding Window
    %           for teeth Teeth
    %--------------------------------------------

    % Initialize arrays to collect windows analysis results
    teeth_x_curr = [];
    teeth_y_curr = [];

    teeth_x_teeth_sep = [];
    teeth_y_teeth_sep = [];

    % The "user-assisted initialization" referenced in the paper
    %   will be interpreted to be the lower intensity detected
    upper_int_mean = [];
    lower_int_mean = [];

    % From paper:
    %    pvi_Di_yi = pvi_Di * pvi_yi
    %    pvi_Di = c * {1 - Di / (max_k * Dk)}
    %    pvi_yi = {1 / [sigma*sqrt(2*pi)]} * exp{-(yi-y_int)^2/(sigma^2)}
    teeth_pvi_yi = [];
    teeth_pvi_Di = [];
    teeth_pvi_Di_yi = [];

    % Iterate through image 
    for teeth_x_curr = 4:(x_im-4)

        %------------------------------------------------------------------
        % Assumption / Approach: 
        % 1- We will iterate through each point on the x-axis
        %    Calculate the mean resolution of the perpendicular line at
        %    that point AND
        %    use that perpendicular line intensity to identify valleys
        %------------------------------------------------------------------
        
        %------------------------------------------------------------------
        % Find the matching y point on the curb
        %------------------------------------------------------------------
        teeth_y_curr = [(teeth_x_curr^2) * coeff_gauss_valley(1) + ...
                        (teeth_x_curr^1) * coeff_gauss_valley(2) + ...
                        (teeth_x_curr^0) * coeff_gauss_valley(3)];

        teeth_y_before = [((teeth_x_curr-1)^2)*coeff_gauss_valley(1)+...
                          ((teeth_x_curr-1)^1)*coeff_gauss_valley(2)+...
                          ((teeth_x_curr-1)^0)*coeff_gauss_valley(3)];

        teeth_y_after = [((teeth_x_curr+1)^2)*coeff_gauss_valley(1)+...
                         ((teeth_x_curr+1)^1)*coeff_gauss_valley(2)+...
                         ((teeth_x_curr+1)^0)*coeff_gauss_valley(3)];

        %------------------------------------------------------------------
        % Find perpendicular to the point in the curb
        %------------------------------------------------------------------

        % Windows / Line Direction of the 2nd Deg Approx
        wind_slope = (teeth_y_after - teeth_y_before) / 2 ;
        
        % Tangent to Windows / Line
        wind_tan_slope = -1 / wind_slope;

        % Find the offset b of the tangent line
        % teeth_y_curr = wind_tan_slope * teeth_x_curr + b
        b = teeth_y_curr - wind_tan_slope * teeth_x_curr;
        %------------------------------------------------------------------


        %------------------------------------------------------------------
        % Find upper and lower bounds perpendicular to the current point
        %------------------------------------------------------------------

        % Upper Limit
        % teeth_y_curr = wind_tan_slope * teeth_x_curr + b
        upper_lim_y = 1;
        upper_lim_x = round((upper_lim_y - b) / wind_tan_slope);

        % If limits outside the picture boundaries
        if (1 > upper_lim_x)
            upper_lim_x =  1;
            upper_lim_y = round(wind_tan_slope * upper_lim_x + b);
        elseif (upper_lim_x > x_im)
            upper_lim_x =  x_im;
            upper_lim_y = round(wind_tan_slope * upper_lim_x + b);
        end

        % Lower Limit
        % teeth_y_curr = wind_tan_slope * teeth_x_curr + b
        lower_lim_y = y_im;
        lower_lim_x = round((lower_lim_y - b) / wind_tan_slope);

        % If limits outside the picture boundaries
        if (1 > lower_lim_x)
            lower_lim_x =  1;
            lower_lim_y = round(wind_tan_slope * lower_lim_x + b);
        elseif (lower_lim_x > x_im)
            lower_lim_x =  x_im;
            lower_lim_y = round(wind_tan_slope * lower_lim_x + b);
        end

        % Save Limits (for plotting)
        % Format to be used: [upper, current, lower]
        teeth_x_teeth_sep = [teeth_x_teeth_sep; [upper_lim_x,...
                                                 teeth_x_curr,...
                                                 lower_lim_x]];

        teeth_y_teeth_sep = [teeth_y_teeth_sep; [upper_lim_y,...
                                                 teeth_y_curr,...
                                                 lower_lim_y]];
        %------------------------------------------------------------------

        %------------------------------------------------------------------
        % Calculate Sweep Line
        %------------------------------------------------------------------
        % Assumption / Approach: 
        % Sweep will have decimal numbers that will be converted to int
        %   this calculation will create the evaluation of duplicate points
        % Assuming there are more duplicated of a specific pixel, this is due
        %   to the fact that more "decimals" are in that pixel, and the
        %   average intensity will therefore give more weight to those values

        % Upper Sweep Values 
        upper_teeth_x_sweep = linspace(teeth_x_curr, ...
                                       upper_lim_x,...
                                       horiz_sample_pts);
        upper_teeth_y_sweep = (wind_tan_slope * upper_teeth_x_sweep + b);


        % Lower Sweep Values 
        lower_teeth_x_sweep = linspace(teeth_x_curr, ...
                                       lower_lim_x,...
                                       horiz_sample_pts);
        lower_teeth_y_sweep = (wind_tan_slope * lower_teeth_x_sweep + b);
        %------------------------------------------------------------------

        %------------------------------------------------------------------
        % Pre-image data massaging
        % Let's trim values that are out of bounds
        %------------------------------------------------------------------

        % Round value so we can use the to get actual pixel intensity
        upper_teeth_x_sweep = round(upper_teeth_x_sweep); 
        upper_teeth_y_sweep = round(upper_teeth_y_sweep); 
        lower_teeth_x_sweep = round(lower_teeth_x_sweep); 
        lower_teeth_y_sweep = round(lower_teeth_y_sweep); 

        % Define sweep and trim out of image values       
        upper_teeth_x_sweep = upper_teeth_x_sweep(x_im > upper_teeth_x_sweep); 
        upper_teeth_y_sweep = upper_teeth_y_sweep(x_im > upper_teeth_x_sweep); 
        lower_teeth_x_sweep = lower_teeth_x_sweep(x_im > lower_teeth_x_sweep); 
        lower_teeth_y_sweep = lower_teeth_y_sweep(x_im > lower_teeth_x_sweep); 

        upper_teeth_x_sweep = upper_teeth_x_sweep(upper_teeth_x_sweep > 0); 
        upper_teeth_y_sweep = upper_teeth_y_sweep(upper_teeth_x_sweep > 0); 
        lower_teeth_x_sweep = lower_teeth_x_sweep(lower_teeth_x_sweep > 0); 
        lower_teeth_y_sweep = lower_teeth_y_sweep(lower_teeth_x_sweep > 0); 

        upper_teeth_x_sweep = upper_teeth_x_sweep(y_im > upper_teeth_y_sweep); 
        upper_teeth_y_sweep = upper_teeth_y_sweep(y_im > upper_teeth_y_sweep); 
        lower_teeth_x_sweep = lower_teeth_x_sweep(y_im > lower_teeth_y_sweep); 
        lower_teeth_y_sweep = lower_teeth_y_sweep(y_im > lower_teeth_y_sweep); 

        upper_teeth_x_sweep = upper_teeth_x_sweep(upper_teeth_y_sweep > 0); 
        upper_teeth_y_sweep = upper_teeth_y_sweep(upper_teeth_y_sweep > 0); 
        lower_teeth_x_sweep = lower_teeth_x_sweep(lower_teeth_y_sweep > 0); 
        lower_teeth_y_sweep = lower_teeth_y_sweep(lower_teeth_y_sweep > 0); 

        %------------------------------------------------------------------

        %------------------------------------------------------------------
        % Let's collect the individual teeth intensity        
        %------------------------------------------------------------------
        
        % Upper teeth - Collect mean intensity
        for upper_tooth = 1:size(upper_teeth_x_sweep,2)
            upper_y_int_mean = mean(...
                                   Image_A(upper_teeth_y_sweep(upper_tooth),...
                                           upper_teeth_x_sweep(upper_tooth))...
                                    );
        end

        % Lower teeth - Collect mean intensity
        for lower_tooth = 1:size(lower_teeth_x_sweep,2)
            lower_y_int_mean = mean(...
                                   Image_A(lower_teeth_y_sweep(lower_tooth),...
                                           lower_teeth_x_sweep(lower_tooth))...
                                    );
        end
            
        % Smooth Intensity
        upper_y_int_mean = smoothn(upper_y_int_mean);
        lower_y_int_mean = smoothn(lower_y_int_mean);

        % Store Data
        upper_int_mean = [upper_int_mean; upper_y_int_mean];
        lower_int_mean = [lower_int_mean; lower_y_int_mean];

    end % End of teeth While Loop

    %--------------------------------------------
    % Step 7.1: Teeth Separation & Plotting
    %--------------------------------------------

    figure(2)
    plot(upper_int_mean,'g-','LineWidth',2)

    figure(3)
    plot(lower_int_mean,'m-','LineWidth',2)

    % Find the local minima
    upper_min_idx = islocalmin(upper_int_mean,...
                               'MinSeparation',55,...
                               'MinProminence',numb_of_teeth);
    lower_min_idx = islocalmin(lower_int_mean,...
                               'MinSeparation',55,...
                               'MinProminence',numb_of_teeth);

    % Upper - Massage Data for plotting
    upper_int_mean_plot = upper_int_mean(upper_min_idx);
    upper_x = teeth_x_teeth_sep(upper_min_idx,[1 2]);
    upper_y = teeth_y_teeth_sep(upper_min_idx,[1 2]);

    % Plot Image
    for upper_tooth = 1:size(upper_x,1)
        figure(1)
        hold on
        plot(upper_x(upper_tooth,:),...
             upper_y(upper_tooth,:),...
             'g-','LineWidth',2)
    end

    % Upper - Massage Data for plotting
    lower_int_mean_plot = lower_int_mean(lower_min_idx);
    lower_x = teeth_x_teeth_sep(lower_min_idx,[2 3]);
    lower_y = teeth_y_teeth_sep(lower_min_idx,[2 3]);

    % Plot Image
    for lower_tooth = 1:size(lower_x,1)
        figure(1)
        hold on
        plot(lower_x(lower_tooth,:),...
             lower_y(lower_tooth,:),...
             'c-','LineWidth',2)
    end

    figure(2)
    plot(upper_int_mean,'g-','LineWidth',2)
    hold on
    plot(teeth_x_teeth_sep(upper_min_idx,2),...
         upper_int_mean(upper_min_idx),...
         'ro','LineWidth',2)
    grid on
    title("Upper Teeth Mean Intensity")

    figure(3)
    plot(lower_int_mean,'c-','LineWidth',2)
    hold on
    plot(teeth_x_teeth_sep(lower_min_idx,2),...
         lower_int_mean(lower_min_idx),...
         'ro','LineWidth',2)
    grid on
    title("Loweer Teeth Mean Intensity")

end         % All images Iteration