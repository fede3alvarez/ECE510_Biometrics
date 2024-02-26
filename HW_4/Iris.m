close all
clear all;
clc;

% Code will be available at
% https://github.com/fede3alvarez/ECE510_Biometrics/HW_4
% Homework 4 - Iris

Available_images = [...
                    'iris1.png'
                    'iris2.png'
                    'iris3.png'
                    'iris4.png'
                    'iris5.png'
                    'iris6.png'
                    ];

f_figure = 1;
gauss = 3;

theta = 0 : 0.01 : 2*pi;

for m = 1:size(Available_images,1)
    current_image = Available_images(m,:);
    [Image_A_Raw, Map_A] = imread(current_image);

    % Filter Image and get Binary
    Image_A = imgaussfilt(Image_A_Raw,gauss);  
        
    % Note: 
    % This script assumes that the Hough_Transform.m script was run first
    % and therefore, pseudo-accurate eye center information is ready 

    % Max and Min Ranges to scan / calculate for pupil and iris    
    min_pupil = 20;
    max_pupil = 50;

    % Load saved data
    load_selections = strcat(current_image(1:5),'_data.mat');
    selections = importdata(load_selections);
    x_center = selections.x_center;
    y_center = selections.y_center;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %               Data  Processing                %
    %                Pupil Analisys                 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Initialize empty matrix to store data
    %       Col 1 = Radius
    %       Col 2 = Ave Intensity
    %       Col 3 = Ave Intensity Gradient
    %       Col 4 = Ave Intensity Hessian
    radius_pupil_data = zeros(floor(max_pupil),4);


    %----------------------------------------------    
    % Iterate through every pupil radius and
    %   get the average intensity for each radius
    %----------------------------------------------
    for radius = 1:max_pupil
    
        % Calculate pupil circle
        x_circle = round(radius * cos(theta) + x_center)';
        y_circle = round(radius * sin(theta) + y_center)';
        samples_idx = [x_circle, y_circle];
        samples_idx = unique(samples_idx, 'rows', 'stable');

        % Get intensity sample
        samples_value = zeros(size(samples_idx,1),1);

        for sample = 1:size(samples_idx,1)

            if (size(Image_A,1) >= samples_idx(sample,1)) &...
               (samples_idx(sample,1) >= 1) &...
               (size(Image_A,2) >= samples_idx(sample,2)) & ...
               (samples_idx(sample,2) >= 1)

                samples_value(sample) = Image_A(...
                                             samples_idx(sample,1),...
                                             samples_idx(sample,2));
            else
                samples_value(sample) = 0;
            end

       end

        % Smoothn values before analysis
        %samples_value = smoothn(samples_value);
        
        % Calculate mean intensity at that radius
        ave_int_pupil = median(samples_value);

        % Store data
        radius_pupil_data(radius,[1 2]) = [radius ave_int_pupil];
    
    end % End of Radius Iterations


    %----------------------------------------------    
    % From the radius average intensity,
    % calculate the gradient, 
    % and find the edge of the pupil
    %----------------------------------------------

    % Calculate the gradient or Derivation
    % for delta = 2:size(radius_pupil_data,1)
    %     radius_pupil_data(delta,3) = radius_pupil_data(delta,2) -...
    %                                 radius_pupil_data(delta-1,2)
    % end
    radius_pupil_data(:,3) = gradient(radius_pupil_data(:,2));
    radius_pupil_data(:,3) = smoothn(radius_pupil_data(:,3));

    % Calculate the Hessian or 2nd derviative
    % for delta = 2:size(radius_pupil_data,1)
    %     radius_pupil_data(delta,4) = radius_pupil_data(delta,3) -...
    %                                 radius_pupil_data(delta-1,3)
    % end

    radius_pupil_data(:,4) = gradient(radius_pupil_data(:,3));
    radius_pupil_data(:,4) = smoothn(radius_pupil_data(:,4));

    % Remove the "Min Radius" part of the data
    radius_pupil_data = radius_pupil_data([min_pupil:max_pupil],:);

    % Select Radius
    [max_grad_pupil, max_grad_pupil_idx] = ...
                     max(radius_pupil_data(:,3));

    %----------------------------------------------    
    % Plop all the graph
    % (Intensity, Gradient and Hessian)
    %----------------------------------------------
    
    figure(f_figure)
    subplot(4,2,5)
    hold on 
    plot(radius_pupil_data(:,1),radius_pupil_data(:,2),'*b-')
    title_pupil_int = strcat(current_image(1:5),', Pupil Intensity');
    title(title_pupil_int)
    hold on

    figure(f_figure)
    subplot(4,2,7)
    hold on 
    plot(radius_pupil_data(:,1),radius_pupil_data(:,3),'*r-')
    title_pupil_grad = strcat(current_image(1:5),', Pupil Gradient');
    title(title_pupil_grad)
    hold on

    % figure(f_figure)
    % subplot(4,2,5)
    % hold on 
    % plot(radius_pupil_data(:,1),radius_pupil_data(:,4),'*g-')
    % title_pupil_grad = strcat(current_image(1:5),', Pupil Hessian');
    % title(title_pupil_grad)
    % hold on

    % Collect data
    [max_pupil_value, max_pupil_idx] = max(radius_pupil_data(:,3));
    pupil_radius = radius_pupil_data(max_pupil_idx,1);
    pupil_int = radius_pupil_data(max_pupil_idx,2);
    pupil_int_grad = radius_pupil_data(max_pupil_idx,3);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %               Data  Processing                %
    %                Iris Analisys                 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Guess ranges
    min_iris = round(pupil_radius*2);
    max_iris = round(pupil_radius*4);

    % Initialize empty matrix to store data
    %       Col 1 = Radius
    %       Col 2 = Ave Intensity
    %       Col 3 = Ave Intensity Gradient
    %       Col 4 = Ave Intensity Hessian
    radius_iris_data = zeros(size(min_iris:1:max_iris,2),4);
   
    
    %----------------------------------------------    
    % Iterate through every pupil radius and
    %   get the average intensity for each radius
    %----------------------------------------------
    for radius = 1:max_iris
    
        % Calculate pupil circle
        x_circle = round(radius * cos(theta) + x_center)';
        y_circle = round(radius * sin(theta) + y_center)';
        samples_idx = [x_circle, y_circle];
        samples_idx = unique(samples_idx, 'rows', 'stable');

        % Get intensity sample
        samples_value = zeros(size(samples_idx,1),1);

        for sample = 1:size(samples_idx,1)

            if (size(Image_A,1) >= samples_idx(sample,1)) &...
               (samples_idx(sample,1) >= 1) &...
               (size(Image_A,2) >= samples_idx(sample,2)) & ...
               (samples_idx(sample,2) >= 1)

                samples_value(sample) = Image_A(...
                                             samples_idx(sample,1),...
                                             samples_idx(sample,2));
            else
                samples_value(sample) = 0;
            end

       end

        % Smoothn values before analysis
        samples_value = smoothn(samples_value);
        
        % Calculate mean intensity at that radius
        ave_int_iris = mean(samples_value);

        % Store data
        radius_iris_data(radius,[1 2]) = [radius ave_int_iris];
    
    end % End of Radius Iterations


    %----------------------------------------------    
    % From the radius average intensity,
    % calculate the gradient, 
    % and find the edge of the pupil
    %----------------------------------------------

    % Calculate the gradient or Derivation
    % for delta = 2:size(radius_iris_data,1)
    %     radius_iris_data(delta,3) = radius_iris_data(delta,2) -...
    %                                 radius_iris_data(delta-1,2)
    % end
    radius_iris_data(:,3) = gradient(radius_iris_data(:,2));
    radius_iris_data(:,3) = smoothn(radius_iris_data(:,3));

    
    % Calculate the Hessian or 2nd derviative
    % for delta = 2:size(radius_iris_data,1)
    %     radius_iris_data(delta,3) = radius_iris_data(delta,3) -...
    %                                 radius_iris_data(delta-1,3)
    % end
    radius_iris_data(:,4) = gradient(radius_iris_data(:,3));
    radius_iris_data(:,4) = smoothn(radius_iris_data(:,4));


    % Select Radius and disregard pupil radius area
    [max_grad_iris, max_grad_iris_idx] = ...
                     max(radius_iris_data([min_iris:end],3));

    % Correct / account for disregarded pupil radius area
    max_grad_iris_idx = max_grad_iris_idx + min_iris;

    iris_radius = radius_iris_data(max_grad_iris_idx,1);
    iris_int = radius_iris_data(max_grad_iris_idx,2);
    iris_int_grad = radius_iris_data(max_grad_iris_idx,3);

    %----------------------------------------------    
    % Plop all the graph
    % (Intensity, Gradient and Hessian)
    %----------------------------------------------

    figure(f_figure)
    subplot(4,2,6)
    hold on 
    plot(radius_iris_data(:,1),radius_iris_data(:,2),'*b-')
    title_iris_int = strcat(current_image(1:5),', Iris Intensity');
    title(title_iris_int)
    hold on

    figure(f_figure)
    subplot(4,2,8)
    hold on 
    plot(radius_iris_data(:,1),radius_iris_data(:,3),'*r-')
    title_iris_grad = strcat(current_image(1:5),', Iris Gradient');
    title(title_iris_grad)
    hold on

    % figure(f_figure)
    % subplot(4,2,6)
    % hold on 
    % plot(radius_iris_data(:,1),radius_iris_data(:,4),'*g-')
    % title_iris_grad = strcat(current_image(1:5),', Iris Hessian');
    % title(title_iris_grad)
    % hold on


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Calculate pupil circle
    x_pupil_plot = pupil_radius * cos(theta) + x_center;
    y_pupil_plot = pupil_radius * sin(theta) + y_center;

    % Calculate iris circle
    x_iris_plot = iris_radius * cos(theta) + x_center;
    y_iris_plot = iris_radius * sin(theta) + y_center;

    % Plot Image
    figure(f_figure)
    subplot(4,2,[1:2])
    imshow(Image_A_Raw);
    title(current_image);    

    % Plot Pupil
    hold on
    plot(x_pupil_plot,y_pupil_plot,'-m');
    hold on

    % Plot Pupil
    hold on
    plot(x_iris_plot,y_iris_plot,'-m');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %            Processing, Save, &                %
    %            Display Iris Matrix                %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Setup the Matrix
    iris_Matrix = zeros((iris_radius - pupil_radius), size(theta,2));

    for angle = 1:size(theta,2)
        ang = theta(angle); 
        for distance = 1:(iris_radius - pupil_radius)
            
            dist = pupil_radius + distance; 
            x_iris_Matrix = round(dist * cos(ang) + x_center);
            y_iris_Matrix = round(dist * sin(ang) + y_center);

            if (size(Image_A_Raw,1) >= x_iris_Matrix) & ...
               (x_iris_Matrix >= 1) & ...
               (size(Image_A_Raw,2) >= y_iris_Matrix) & ...
               (y_iris_Matrix >= 1)
                iris_Matrix(distance,angle) = ...
                    Image_A_Raw(x_iris_Matrix, y_iris_Matrix);
            end % End-if to populate Matrix
       end % End of Radius Loop
    end % End of Angle Loop

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Plot Image
    subplot(4,2,[3:4])
    figure(f_figure)
    imshow(iris_Matrix);
    iris_title = strcat('Iris Matrix-',current_image);
    title(iris_title);
    f_figure = f_figure + 1;

end % End of Images Iterations