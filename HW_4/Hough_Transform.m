close all
clear all;
clc;

% Code will be available at
% https://github.com/fede3alvarez/ECE510_Biometrics/HW_4
% Homework 4 - Hough Transform

Available_images = [...
                    'iris1.png'
                    'iris2.png'
                    'iris3.png'
                    'iris4.png'
                    'iris5.png'
                    'iris6.png'
                    ];

f_figure = 1;
gauss = 5;

for m = 1:size(Available_images,1)
    current_image = Available_images(m,:);
    [Image_A, Map_A] = imread(current_image);

    % Filter Image and get Binary
    Image_A = imgaussfilt(Image_A,gauss);  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                Finding the              %
    %             Gradients Crossing          %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Initialize Matrix to keep track gradient x-ing
    Grad_Crossroads = zeros(size(Image_A));
    [Grad_mag, Grad_dir] = imgradient(Image_A);
    Grad_dir = deg2rad(Grad_dir);
    
    % Helper Values  
    x_steps = zeros(size(Image_A,1),1);
    y_steps = zeros(size(Image_A,1),1);
    % Old fashion iteration per element
    for row = 1:size(Image_A,1)
        for col = 1:size(Image_A,2)

            % Before doing any math, 
            %   less make sure the gradient is not zero
            if (Grad_mag(row,col) ~= 0)
            
                % Helper Values - to keep track of repeated line entries
                x_step_old = -inf;
                y_step_old = -inf;

                % If Gradient is not zero, 
                %   go through each x value in image
                %   and calculate the gradient line / direction
                for x_step = 1:size(Image_A,1);
                    % Calculate match y-value
                    y_step = round((x_step - row)* tan(Grad_dir(row,col)+pi/2) + col);

                    % Populate Map if the value is not a repeat
                    %   and value is within bounds
                    if ((x_step_old ~= x_step) &...
                        (y_step_old ~= y_step) &...
                        (size(Image_A,1) >= x_step) &...
                        (size(Image_A,2) >= y_step) &...
                        (x_step >= 1) & (y_step >= 1))
                            Grad_Crossroads(x_step,y_step) = ...
                                Grad_Crossroads(x_step,y_step) + 1;
                    end % End of updating Map Value

                    % Update Helper Values
                    x_step_old = x_step;
                    y_step_old = y_step;
                    
                    % Update Helper Values
                    x_steps(x_step,1) = x_step;
                    y_steps(x_step,1) = y_step;

                end % End- Line Calculation

            end % End-If
        end % End-Col
    end % End-Row


    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %              Plot the Center            %
    %                 Candidates              %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Find Center candidates based on the max accumulation
    [y_center, x_center] = find(ismember(Grad_Crossroads,...
                                         max(Grad_Crossroads(:))));

 
    % Save data
    save_file_selections = strcat(current_image(1:5),'_data')
    save(save_file_selections,'x_center','y_center')

    % Plot Center Candidates
    pause(3);
    figure(f_figure)
    imshow(Image_A);
    title("Center Candidates");
    hold on
    plot(x_center,y_center,'*m','LineWidth',2);
    f_figure = f_figure + 1;
    pause(3);

end % End-Image