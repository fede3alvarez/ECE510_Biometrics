close all
clear all;
clc;

% Code will be available at
% https://github.com/fede3alvarez/ECE510_Biometrics/HW_4
% Homework 4 - Iris

Available_images = ['iris1.jpg'
                    % 'iris2.jpg'
                    % 'iris3.jpg'
                    % 'iris4.jpg'
                    % 'iris5.jpg'
                    % 'iris6.jpg'
                    ];

f_figure = 1;

% For each point find the Feature Metric / Finger Thickness
feature_metric = [];

for m = 1:size(Available_images,1)
    current_image = Available_images(m,:);
    [Image_A, Map_A] = imread(current_image);
    
    % Filter Image
    Image_A = imgaussfilt(Image_A,7);
    
    % % Plot Image
    % figure(f_figure)
    % %subplot(1,1,f_figure)
    % imshow(Image_A);
    % title(current_image);
        
    % Note: 
    % The user selected an approximated eyes center 
    % before running this script

    % Max and Min Ranges to scan / calcualte for pupil and iris
    % Note: Ratios were obtained by trial and error on the given images

    % No min_pupil is needed: this is the center of the eye
    max_pupil = size(Image_A,1)/8;

    min_iris = size(Image_A,1)/16;
    max_iris = size(Image_A,1)/2;

    % Load saved data
    load_selections = strcat(current_image(1:5),'_data.mat');
    selections = importdata(load_selections);
    x = selections.x;
    y = selections.y;

    x_raw_center = x;
    y_raw_center = y;

    % % Plot Eye Center
    % hold on
    % plot(x,y,'*m','LineWidth',3);

    % % Debuging plots
    % % Plot Pupil Max Range
    % hold on
    % plot([x (x+max_pupil)],[y y],'*-c');
    % 
    % % Plot Iris Max and Min Ranges
    % hold on
    % plot([(x-max_iris) (x-min_iris)],[y y],'*-g');
    % 
    % center_radius = 15;
    % center_box = [(x - center_radius) (y + center_radius)
    %               (x + center_radius) (y + center_radius)
    %               (x + center_radius) (y - center_radius)
    %               (x - center_radius) (y - center_radius)
    %               (x - center_radius) (y + center_radius)];
    % 
    % % Plot Eye Center Boxes
    % hold on
    % plot(center_box(:,1), center_box(:,2),'*-g');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %               Data  Processing                %
    %                Pupil Analisys                 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % From the approximate raw center specified by the user
    % Create a "box" of center_radius of possible center candidates
    center_radius = 3;   
    [x_center_box, y_center_box] = meshgrid(...
                              [(floor(x_raw_center) - center_radius):...
                               1:...
                               (ceil(x) + center_radius)],...
                              [(floor(y_raw_center) - center_radius):...
                               1:...
                               (ceil(y_raw_center) + center_radius)]);

    center_box = [reshape(x_center_box,...
                          size(x_center_box,1)*size(x_center_box,2),...
                          1),...
                  reshape(y_center_box,...
                          size(y_center_box,1)*size(y_center_box,2),...
                          1)];

    % Initialize empty matrix to store data
    center_data = zeros(size(center_box,1),5);

    % Iterate through each center candidate
    % The center with the largest pupil gradient will be 
    % consider the final center
    for eye_center = 1:size(center_box,1)
        x_eye_center = center_box(eye_center,1);
        y_eye_center = center_box(eye_center,2);
    
        % Initialize empty matrix to store data
        radius_pupil_data = zeros(floor(max_pupil),3);
        radius_pupil_thickness = 0;
    
            % Iterate through every pupil radius
            for radius = 1:max_pupil
            
                % Create a circle map image to identify the pupil circle
                [pupil_cols,pupil_rows] = meshgrid(1:size(Image_A,1),...
                                                   1:size(Image_A,2));
                circlePixels = ((pupil_rows - x_eye_center).^2 +...
                                (pupil_cols - y_eye_center).^2 ...
                             <= (radius.^2 + radius_pupil_thickness)) & ...
                               ((pupil_rows - x_eye_center).^2 +...
                                (pupil_cols - y_eye_center).^2 ...
                             >= (radius.^2 - radius_pupil_thickness));
            
                % Get intensity at circle
                circlePixels = Image_A(circlePixels);
            
                % Calculate average intensity at that radius
                ave_int_pupil = mean(circlePixels);
            
                % Store data
                radius_pupil_data(radius,[1 2]) = [radius ave_int_pupil];
            end % End of Radius Iterations
    
        radius_pupil_data(:,3) = gradient(radius_pupil_data(:,2));
        [max_grad_pupil, max_grad_pupil_idx] = max(radius_pupil_data(:,3));

        center_data(eye_center,:) = [x_eye_center,...
                                 y_eye_center,...
                                 radius_pupil_data(max_grad_pupil_idx,:)];
    
        figure(2)
        hold on 
        plot(radius_pupil_data(:,1),radius_pupil_data(:,2),'b-')
        title("Intensity")
        hold on

        figure(3)
        hold on 
        plot(radius_pupil_data(:,1),radius_pupil_data(:,3),'r-')
        title("Gradient")
        hold on

    end % End of Centers Iterations

    % Collect data for final center and pupil radius
    [max_pupil_value, max_pupil_idx] = max(center_data(:,5));
    x_final_center = center_data(max_pupil_idx,1);
    y_final_center = center_data(max_pupil_idx,2);
    pupil_radius = center_data(max_pupil_idx,3);
    pupil_int = center_data(max_pupil_idx,4);
    pupil_int_grad = center_data(max_pupil_idx,5);

    % Plot Image
    figure(1)
    %subplot(1,1,f_figure)
    imshow(Image_A);
    title(current_image);

    % Plot Eye Center
    figure(1)
    hold on
    plot(x_final_center,y_final_center,'*m','LineWidth',2);

    % Calculate pupil circle
    theta = 0 : 0.01 : 2*pi;
    x_pupil_plot = pupil_radius * cos(theta) + x_final_center;
    y_pupil_plot = pupil_radius * sin(theta) + y_final_center;

    % Plot Pupil
    figure(1)
    hold on
    plot(x_pupil_plot,y_pupil_plot,'-m');

   %end % End of Centers Iterations 

end % End of Images Iterations