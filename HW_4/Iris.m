close all
clear all;
clc;

% Code will be available at
% https://github.com/fede3alvarez/ECE510_Biometrics/HW_4
% Homework 4 - Iris

Available_images = [...
                    %'iris1.jpg'
                     'iris2.jpg'
                    % 'iris3.jpg'
                    % 'iris4.jpg'
                    % 'iris5.jpg'
                    % 'iris6.jpg'
                    ];

f_figure = 1;
gauss = 3;

theta = 0 : 0.01 : 2*pi;

for m = 1:size(Available_images,1)
    current_image = Available_images(m,:);
    [Image_A, Map_A] = imread(current_image);
    
    if not(ismatrix(Image_A))
        Image_A = rgb2gray(Image_A);
    end

    % Filter Image
    %Image_A = imgaussfilt(Image_A,gauss);
        
    % Note: 
    % The user selected an approximated eyes center 
    % before running this script

    % Max and Min Ranges to scan / calcualte for pupil and iris
    
    min_pupil = 20;
    max_pupil = size(Image_A,1)/7;

    % Load saved data
    load_selections = strcat(current_image(1:5),'_data.mat');
    selections = importdata(load_selections);
    x = selections.x;
    y = selections.y;

    x_raw_center = x;
    y_raw_center = y;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %               Data  Processing                %
    %                Pupil Analisys                 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % From the approximate raw center specified by the user
    % Create a "box" of center_radius of possible center candidates
    center_radius = 0;   
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
    center_data = zeros(size(center_box,1),6);
    

    % Iterate through each center candidate
    % The center with the largest pupil gradient will be 
    % consider the final center
    for eye_center = 1:size(center_box,1)
        x_eye_center = center_box(eye_center,1);
        y_eye_center = center_box(eye_center,2);
    
        % Initialize empty matrix to store data
        radius_pupil_data = zeros(floor(max_pupil),4);
        
            % Iterate through every pupil radius
            for radius = 1:max_pupil
            
                % Calculate pupil circle
                theta = 0 : 0.01 : 2*pi;
                x_circle = round(radius * cos(theta) + x_eye_center)';
                y_circle = round(radius * sin(theta) + y_eye_center)';
                samples_idx = [x_circle, y_circle];
                samples_idx = unique(samples_idx, 'rows', 'stable');

                % Get intensity sample
                samples_value = zeros(size(samples_idx,1),1);

                for sample = 1:size(samples_idx,1)
                    samples_value(sample) = Image_A(...
                                                 samples_idx(sample,1),...
                                                 samples_idx(sample,2));
                end

                % Smoothn values before analysis
                samples_value = smoothn(samples_value);
                
                % Calculate mean intensity at that radius
                ave_int_pupil = mean(samples_value);

                % Store data
                radius_pupil_data(radius,[1 2]) = [radius ave_int_pupil];
            
            end % End of Radius Iterations
    
        % Calculate the gradient or Derivation
        for delta = 2:size(radius_pupil_data,1)
            radius_pupil_data(delta,3) = radius_pupil_data(delta,2) -...
                                        radius_pupil_data(delta-1,2)
        end
        radius_pupil_data(:,3) = smoothn(radius_pupil_data(:,3));

        % Calculate the Hessian or 2nd derviative
        for delta = 2:size(radius_pupil_data,1)
            radius_pupil_data(delta,4) = radius_pupil_data(delta,3) -...
                                        radius_pupil_data(delta-1,3)
        end
        radius_pupil_data(:,4) = smoothn(radius_pupil_data(:,4));

        % Remove the "Min Radius" part of the data
        radius_pupil_data = radius_pupil_data([min_pupil:max_pupil],:);

        % Select Radius
        [max_grad_pupil, max_grad_pupil_idx] = ...
                         max(radius_pupil_data(:,3));

        % Store Data
        center_data(eye_center,:) = [x_eye_center,...
                                 y_eye_center,...
                                 radius_pupil_data(max_grad_pupil_idx,:)];
        
        figure(f_figure)
        subplot(3,2,1)
        hold on 
        plot(radius_pupil_data(:,1),radius_pupil_data(:,2),'*b-')
        title_pupil_int = strcat(current_image(1:5),', Pupil Intensity');
        title(title_pupil_int)
        hold on

        figure(f_figure)
        subplot(3,2,3)
        hold on 
        plot(radius_pupil_data(:,1),radius_pupil_data(:,3),'*r-')
        title_pupil_grad = strcat(current_image(1:5),', Pupil Gradient');
        title(title_pupil_grad)
        hold on

        figure(f_figure)
        subplot(3,2,5)
        hold on 
        plot(radius_pupil_data(:,1),radius_pupil_data(:,4),'*g-')
        title_pupil_grad = strcat(current_image(1:5),', Pupil Hessian');
        title(title_pupil_grad)
        hold on


    end % End of Centers Iterations

    % Collect data for final center and pupil radius
    [max_pupil_value, max_pupil_idx] = max(center_data(:,5));
    x_final_center = center_data(max_pupil_idx,1);
    y_final_center = center_data(max_pupil_idx,2);
    pupil_radius = center_data(max_pupil_idx,3);
    pupil_int = center_data(max_pupil_idx,4);
    pupil_int_grad = center_data(max_pupil_idx,5);



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %               Data  Processing                %
    %                Iris Analisys                 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    min_iris = pupil_radius;
    max_iris = pupil_radius*2.5;

    % Initialize empty matrix to store data
    radius_iris_data = zeros(size(pupil_radius:1:max_iris,2),4);
    
    % Iterate through every pupil radius
    for radius = pupil_radius:1:max_iris
    
        % Calculate pupil circle
        theta = 0 : 0.01 : 2*pi;
        x_circle = round(radius * cos(theta) + x_final_center)';
        y_circle = round(radius * sin(theta) + y_final_center)';
        samples_idx = [x_circle, y_circle];
        samples_idx = unique(samples_idx, 'rows', 'stable');

        % Get intensity sample
        samples_value = zeros(size(samples_idx,1),1);

        for sample = 1:size(samples_idx,1)

            if (samples_idx(sample,1)<= size(Image_A,1)) && ...
               (samples_idx(sample,2)<= size(Image_A,2))
            samples_value(sample) = Image_A(...
                                         samples_idx(sample,1),...
                                         samples_idx(sample,2));
        
            end
        end
        % Smoothn values before analysis
        samples_value = smoothn(samples_value);
        
        % Calculate mean intensity at that radius
        ave_int_iris = mean(samples_value);

        % Store data
        radius_iris_data(radius,[1 2]) = [radius ave_int_iris];
    
    end % End of Radius Iterations

    % Calculate the gradient or Derivation
    for delta = 2:size(radius_iris_data,1)
        radius_iris_data(delta,3) = radius_iris_data(delta,2) -...
                                    radius_iris_data(delta-1,2)
    end
    radius_iris_data(:,3) = smoothn(radius_iris_data(:,3));

    % Calculate the Hessian or 2nd derviative
    for delta = 2:size(radius_iris_data,1)
        radius_iris_data(delta,3) = radius_iris_data(delta,3) -...
                                    radius_iris_data(delta-1,3)
    end
    radius_iris_data(:,4) = smoothn(radius_iris_data(:,4));

    % Remove the "Min Radius" part of the data
    radius_iris_data = radius_iris_data([min_iris:max_iris],:);

    % Select Radius
    [max_grad_iris, max_grad_iris_idx] = ...
                     max(radius_iris_data(:,3));


    iris_radius = radius_iris_data(max_grad_iris_idx,1);
    iris_int = radius_iris_data(max_grad_iris_idx,2);
    iris_int_grad = radius_iris_data(max_grad_iris_idx,3);

    figure(f_figure)
    subplot(3,2,2)
    hold on 
    plot(radius_iris_data(:,1),radius_iris_data(:,2),'*b-')
    title_iris_int = strcat(current_image(1:5),', Iris Intensity');
    title(title_iris_int)
    hold on

    figure(f_figure)
    subplot(3,2,4)
    hold on 
    plot(radius_iris_data(:,1),radius_iris_data(:,3),'*r-')
    title_iris_grad = strcat(current_image(1:5),', Iris Gradient');
    title(title_iris_grad)
    hold on

    figure(f_figure)
    subplot(3,2,6)
    hold on 
    plot(radius_iris_data(:,1),radius_iris_data(:,4),'*g-')
    title_iris_grad = strcat(current_image(1:5),', Iris Hessian');
    title(title_iris_grad)
    hold on


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % Plot Image
    f_figure = f_figure + 1;
    figure(f_figure)
    imshow(Image_A);
    title(current_image);
    
    % Helper Line
    x_help = x_final_center:1:(x_final_center+max_pupil);
    y_help = y_final_center*ones(size(x_help));

    x_help_2 = (x_final_center-min_pupil):1:x_final_center;
    y_help_2 = y_final_center*ones(size(x_help_2));

    % Helper Line
    x_help_3 = x_final_center:1:(x_final_center+max_iris);
    y_help_3 = y_final_center*ones(size(x_help_3));

    x_help_4 = (x_final_center-min_iris):1:x_final_center;
    y_help_4 = y_final_center*ones(size(x_help_4));


    figure(f_figure)
    hold on
    %plot([x_final_center,y_final_center,'*m','LineWidth',2);
    plot(x_help,y_help,'*-b','LineWidth',2);

    figure(f_figure)
    hold on
    %plot([x_final_center,y_final_center,'*m','LineWidth',2);
    plot(x_help_2,y_help_2,'*-r','LineWidth',2);

    figure(f_figure)
    hold on
    %plot([x_final_center,y_final_center,'*m','LineWidth',2);
    plot(x_help_3,y_help_3,'*-g','LineWidth',1);

    figure(f_figure)
    hold on
    %plot([x_final_center,y_final_center,'*m','LineWidth',2);
    plot(x_help_4,y_help_4,'*-c','LineWidth',1);

    % Calculate pupil circle
    theta = 0 : 0.01 : 2*pi;
    x_pupil_plot = pupil_radius * cos(theta) + x_final_center;
    y_pupil_plot = pupil_radius * sin(theta) + y_final_center;

    % Plot Pupil
    figure(f_figure)
    hold on
    plot(x_pupil_plot,y_pupil_plot,'-m');
    hold on

    % Calculate iris circle
    theta = 0 : 0.01 : 2*pi;
    x_iris_plot = iris_radius * cos(theta) + x_final_center;
    y_iris_plot = iris_radius * sin(theta) + y_final_center;

    % Plot Pupil
    figure(f_figure)
    hold on
    plot(x_iris_plot,y_iris_plot,'-m');
    f_figure = f_figure + 1;
    hold on
end % End of Images Iterations