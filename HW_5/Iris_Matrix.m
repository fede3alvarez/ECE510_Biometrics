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

    % Load saved data
    load_selections = strcat(current_image(1:5),'_data.mat');
    selections = importdata(load_selections);
    x_center = selections.x_center;
    y_center = selections.y_center;
    pupil_radius = selections.pupil_radius;
    iris_radius = selections.iris_radius;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %            Processing, Save, &                %
    %            Display Iris Matrix                %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Setup the Matrix
    iris_Matrix = zeros((iris_radius - pupil_radius), size(theta,2));

    for angle = 1:size(theta,2)
         ang = theta(angle); 

         for rad = 1:(iris_radius - pupil_radius)
    
            radius = pupil_radius + rad;
    
            % Calculate pupil circle
            x_circle = round(radius * cos(ang) + x_center)';
            y_circle = round(radius * sin(ang) + y_center)';
            samples_idx = [x_circle, y_circle];
    
            % Iterate through the Iris Circle
            for sample = 1:size(samples_idx,1)
    
                if (size(Image_A_Raw,1) >= samples_idx(sample,1)) &...
                   (samples_idx(sample,1) >= 1) &...
                   (size(Image_A_Raw,2) >= samples_idx(sample,2)) & ...
                   (samples_idx(sample,2) >= 1)
    
                    iris_Matrix(radius,angle) = Image_A_Raw(...
                                                 samples_idx(sample,1),...
                                                 samples_idx(sample,2));
                end
            end
    
               
          end % End of Radius Loop
    end % End of Angle Loop

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Plot Image
    figure(f_figure)
    heatmap(iris_Matrix);
    iris_title = strcat('Iris Matrix-',current_image);
    title(iris_title);
    f_figure = f_figure + 1;

end % End of Images Iterations