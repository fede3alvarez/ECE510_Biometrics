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

c = 1;
gauss = 1;
horiz_sample_pts = 150;
horiz_sample_edge = 0.4;

%---------------------------------------
% Parameters for Individual 
%   Teeth calculation
%---------------------------------------
upper_min_dist = 70;
numb_of_teeth = 5;

lower_min_dist = 15;
lower_numb_of_teeth = 7;

for m = 1:size(Available_images,1)

    %---------------------------------------
    % Step 0: Load and Plot Image
    %---------------------------------------

    % In this case, it is in Grayscale
    current_image = Available_images(m,:);
    [Image_A, Map_A] = imread(current_image);
    % Image_A = imgaussfilt(Image_A,3);  
    % Image_A = imgaussfilt(Image_A,3);  
    
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

    % Plot Image
    figure(1)
    hold on
    plot(x_gauss,y_gauss,'m*')


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

    % Iterate through image 
    for upper_x_curr = 4:(x_im-4)

        % Assumption / Approach: 
        % 1- We will iterate through each point on the x-axis
        %    Calculate the mean resolution of the perpendicular line at
        %    that point AND
        %    use that perpendicular line intensity to identify valleys

        % Find the matching y point on the curb
        upper_y_curr = [(upper_x_curr^2) * coeff_gauss_valley(1) + ...
                        (upper_x_curr^1) * coeff_gauss_valley(2) + ...
                        (upper_x_curr^0) * coeff_gauss_valley(3)];

        upper_y_before = [((upper_x_curr-0.5)^2)*coeff_gauss_valley(1)+...
                          ((upper_x_curr-0.5)^1)*coeff_gauss_valley(2)+...
                          ((upper_x_curr-0.5)^0)*coeff_gauss_valley(3)];

        upper_y_after = [((upper_x_curr+0.5)^2)*coeff_gauss_valley(1)+...
                         ((upper_x_curr+0.5)^1)*coeff_gauss_valley(2)+...
                         ((upper_x_curr+0.5)^0)*coeff_gauss_valley(3)];
        
        % Windows / Line Direction of the 2nd Deg Approx
        wind_slope = (upper_y_after - upper_y_before) / 2 ;
        
        % Tangent to Windows / Line
        wind_tan_slope = -1 / wind_slope;

        % Find the offset b of the tangent line
        % upper_y_curr = wind_tan_slope * upper_x_curr + b
        b = upper_y_curr - wind_tan_slope * upper_x_curr;

        % Find x-intercept
        upper_x_intpt_y = 0;
        upper_x_intpt_x = round((upper_x_intpt_y - b) / wind_tan_slope);
       
        % Find y-intercept
        upper_y_intpt_x = 0;
        upper_y_intpt_y = round(wind_tan_slope * upper_x_intpt_x + b);

        % Calculate Sweep Line
        %    i.e From point in Valley 2nd Degree Approx to where?
        if (upper_y_intpt_x > 0)
            upper_x_lim = upper_y_intpt_x;
            upper_y_lim = upper_y_intpt_y;
        else
            upper_x_lim = upper_x_intpt_x;
            upper_y_lim = upper_x_intpt_y;
        end

        if (upper_x_curr > upper_x_lim)
            upper_x_sweep = round(upper_x_lim+1):round(upper_x_curr-1);
            
            % Store values
            upper_x_teeth_sep = [upper_x_teeth_sep;...
                                     [round(upper_x_lim),...
                                      round(upper_x_curr)]...
                                ];

            upper_y_teeth_sep = [upper_y_teeth_sep;...
                                     [round(upper_y_lim),...
                                      round(upper_y_curr)]...
                                ];

        else
            upper_x_sweep = round(upper_x_curr-1):...
                            round(upper_x_lim+1);

            % Store values
            upper_x_teeth_sep  = [upper_x_teeth_sep;...
                                      [round(upper_x_curr),...
                                       round(upper_x_lim)]...
                                  ];

            upper_y_teeth_sep = [upper_y_teeth_sep;...
                                     [round(upper_y_curr),...
                                      round(upper_y_lim)]...
                                ];
        end

        % Define sweep and trim out of image values
        upper_x_sweep = round(upper_x_sweep(x_im > upper_x_sweep));
        upper_y_sweep = round((upper_x_sweep - b) / wind_tan_slope);
        
        % Ensure sweep range is within image limits
        upper_x_sweep = upper_x_sweep(y_im > upper_y_sweep);
        upper_y_sweep = upper_y_sweep(y_im > upper_y_sweep);



        % Round Edges - NO trimming
        upper_low_lim = round(1); %round(size(upper_x_sweep,2)*0.1)
        upper_high_lim = round(size(upper_x_sweep,2)*1);
        
        if (0 >= upper_low_lim)
            upper_low_lim = 1;
        end

        if (upper_high_lim > size(upper_x_sweep,2))
            upper_high_lim = size(upper_x_sweep,2);
        end

        upper_x_sweep = upper_x_sweep(upper_low_lim:upper_high_lim);
        upper_y_sweep = upper_y_sweep(upper_low_lim:upper_high_lim);

        %unique(sort(upper_sweep,2), 'rows')
        % Sweept values, collect mean, and smooth it
        for upper_sweep_pt = 1:size(upper_x_sweep,2)
            upper_y_int_mean = mean(...
                            Image_A(upper_y_sweep(upper_sweep_pt),...
                                    upper_x_sweep(upper_sweep_pt))...
                                    );
        end
            
        % Store values
        upper_y_int_mean = smoothn(upper_y_int_mean);
        upper_x_int = [upper_x_int; upper_x_curr];
        upper_y_int = [upper_y_int; upper_y_int_mean];


    end % End of Upper While Loop

    %--------------------------------------------
    % Step 5.2: Define Sliding Window
    %           for Lower Teeth
    %--------------------------------------------

      % Initialize arrays to collect windows analysis results
    lower_x_curr = [];
    lower_y_curr = [];

    lower_x_teeth_sep = [];
    lower_y_teeth_sep = [];

    % The "user-assisted initialization" referenced in the paper
    %   will be interpreted to be the lower intensity detected
    lower_x_int = [];
    lower_y_int = [];

    % From paper:
    %    pvi_Di_yi = pvi_Di * pvi_yi
    %    pvi_Di = c * {1 - Di / (max_k * Dk)}
    %    pvi_yi = {1 / [sigma*sqrt(2*pi)]} * exp{-(yi-y_int)^2/(sigma^2)}
    lower_pvi_yi = [];
    lower_pvi_Di = [];
    lower_pvi_Di_yi = [];

    % Iterate through image 
    for lower_x_curr = 4:(x_im-4)

        % Assumption / Approach: 
        % 1- We will iterate through each point on the x-axis
        %    Calculate the mean resolution of the perpendicular line at
        %    that point AND
        %    use that perpendicular line intensity to identify valleys

        % Find the matching y point on the curb
        lower_y_curr = [(lower_x_curr^2) * coeff_gauss_valley(1) + ...
                        (lower_x_curr^1) * coeff_gauss_valley(2) + ...
                        (lower_x_curr^0) * coeff_gauss_valley(3)];

        lower_y_before = [((lower_x_curr-0.5)^2)*coeff_gauss_valley(1)+...
                          ((lower_x_curr-0.5)^1)*coeff_gauss_valley(2)+...
                          ((lower_x_curr-0.5)^0)*coeff_gauss_valley(3)];

        lower_y_after = [((lower_x_curr+0.5)^2)*coeff_gauss_valley(1)+...
                         ((lower_x_curr+0.5)^1)*coeff_gauss_valley(2)+...
                         ((lower_x_curr+0.5)^0)*coeff_gauss_valley(3)];
        
        % Windows / Line Direction of the 2nd Deg Approx
        wind_slope = (lower_y_after - lower_y_before) / 2 ;
        
        % Tangent to Windows / Line
        wind_tan_slope = -1 / wind_slope;

        % Find the offset b of the tangent line
        % lower_y_curr = wind_tan_slope * lower_x_curr + b
        b = lower_y_curr - wind_tan_slope * lower_x_curr;

        % Find x-intercept
        lower_x_intpt_y = y_im;
        lower_x_intpt_x = round((lower_x_intpt_y - b) / wind_tan_slope);
       
        % Find y-intercept
        lower_y_intpt_x = x_im;
        lower_y_intpt_y = round(wind_tan_slope * lower_y_intpt_x + b);

        % Calculate Sweep Line
        %    i.e From point in Valley 2nd Degree Approx to where?
        if (lower_x_intpt_y > x_im)
            lower_x_lim = lower_y_intpt_x;
            lower_y_lim = lower_y_intpt_y;
        else
            lower_x_lim = lower_x_intpt_x;
            lower_y_lim = lower_x_intpt_y;
        end

        if (lower_x_lim > lower_x_curr)
            lower_x_sweep = round(lower_x_lim+1):round(lower_x_curr-1);
            
            % Store values
            lower_x_teeth_sep = [lower_x_teeth_sep;...
                                     [round(lower_x_lim),...
                                      round(lower_x_curr)]...
                                ];

            lower_y_teeth_sep = [lower_y_teeth_sep;...
                                     [round(lower_y_lim),...
                                      round(lower_y_curr)]...
                                ];

        else
            lower_x_sweep = round(lower_x_curr-1):...
                            round(lower_x_lim+1);

            % Store values
            lower_x_teeth_sep  = [lower_x_teeth_sep;...
                                      [round(lower_x_curr),...
                                       round(lower_x_lim)]...
                                  ];

            lower_y_teeth_sep = [lower_y_teeth_sep;...
                                     [round(lower_y_curr),...
                                      round(lower_y_lim)]...
                                ];
        end

        % Define sweep and trim out of image values
        lower_x_sweep = round(lower_x_sweep(x_im > lower_x_sweep));
        lower_y_sweep = round((lower_x_sweep - b) / wind_tan_slope);
        
        % Ensure sweep range is within image limits
        lower_x_sweep = lower_x_sweep(y_im > lower_y_sweep);
        lower_y_sweep = lower_y_sweep(y_im > lower_y_sweep);

        lower_x_sweep = lower_x_sweep(lower_y_sweep > 0);
        lower_y_sweep = lower_y_sweep(lower_y_sweep > 0);


        % Round Edges - NO trimming
        lower_low_lim = round(1); %round(size(lower_x_sweep,2)*0.1)
        lower_high_lim = round(size(lower_x_sweep,2)*1);
        
        if (0 >= lower_low_lim)
            lower_low_lim = 1;
        end

        if (lower_high_lim > size(lower_x_sweep,2))
            lower_high_lim = size(lower_x_sweep,2);
        end

        lower_x_sweep = lower_x_sweep(lower_low_lim:lower_high_lim);
        lower_y_sweep = lower_y_sweep(lower_low_lim:lower_high_lim);

        %unique(sort(lower_sweep,2), 'rows')
        % Sweept values, collect mean, and smooth it
        for lower_sweep_pt = 1:size(lower_x_sweep,2)
            lower_y_int_mean = mean(...
                            Image_A(lower_y_sweep(lower_sweep_pt),...
                                    lower_x_sweep(lower_sweep_pt))...
                                    );
        end
            
        % Store values
        lower_y_int_mean = smoothn(lower_y_int_mean);
        lower_x_int = [lower_x_int; lower_x_curr];
        lower_y_int = [lower_y_int; lower_y_int_mean];


    end % End of lower While Loop
    
    %--------------------------------------------
    % Step 7.1: Teeth Separation
    %--------------------------------------------


    upper_inspection_band = size(upper_y_int,1);

    upper_y_int_mean_plot = upper_y_int;
    upper_x_teeth_sep_plot = upper_x_teeth_sep;
    upper_y_teeth_sep_plot = upper_y_teeth_sep;

    lower_y_int_mean_plot = lower_y_int;
    lower_x_teeth_sep_plot = lower_x_teeth_sep;
    lower_y_teeth_sep_plot = lower_y_teeth_sep;

    % % Plot Image
    f_figure = f_figure + 1;
    figure(f_figure)
    plot(lower_x_int, lower_y_int)

    % Plot Image
    for teeth = 1:numb_of_teeth

        % Find min intensity
        [upper_min_val, upper_min_idx] = min(upper_y_int_mean_plot);
        [lower_min_val, lower_min_idx] = min(lower_y_int_mean_plot);

        % Plot Intensity / Teeth limit
        figure(1)
        hold on
        plot(upper_x_teeth_sep_plot(upper_min_idx,:),...
             upper_y_teeth_sep_plot(upper_min_idx,:),...
             'm-','LineWidth',2)
        hold on
        plot(lower_x_teeth_sep_plot(lower_min_idx,:),...
             lower_y_teeth_sep_plot(lower_min_idx,:),...
             'g-','LineWidth',2)

        upper_clip_start = upper_min_idx - upper_min_dist;
        upper_clip_end = upper_min_idx + upper_min_dist;

        if (1 > upper_clip_start)
            upper_clip_start = 1;
        end

        if (upper_clip_end > size(upper_y_int_mean_plot))
            upper_clip_end = size(upper_y_int_mean_plot);
        end

        % Remove it for next iteration
        upper_y_int_mean_plot(upper_clip_start:upper_clip_end) = [];
        upper_x_teeth_sep_plot(upper_clip_start:upper_clip_end,:) = [];
        upper_y_teeth_sep_plot(upper_clip_start:upper_clip_end,:) = [];

        lower_clip_start = lower_min_idx - lower_min_dist;
        lower_clip_end = lower_min_idx + lower_min_dist;

        if (1 > lower_clip_start)
            lower_clip_start = 1;
        end

        if (lower_clip_end > size(lower_y_int_mean_plot))
            lower_clip_end = size(lower_y_int_mean_plot);
        end

        % Remove it for next iteration
        lower_y_int_mean_plot(lower_clip_start:lower_clip_end) = [];
        lower_x_teeth_sep_plot(lower_clip_start:lower_clip_end,:) = [];
        lower_y_teeth_sep_plot(lower_clip_start:lower_clip_end,:) = [];
    end
end         % All images Iteration