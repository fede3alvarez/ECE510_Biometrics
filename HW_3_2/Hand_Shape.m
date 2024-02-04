close all
clear all;
clc;

% Code will be available at
% https://github.com/fede3alvarez/ECE510_Biometrics
% Homework 3 - Hand Shapes

% Images
Available_images = ['HandShape00.jpg'
                    ];

for m = 1:size(Available_images,1)
    current_image = Available_images(m,:);
    [Image_A, Map_A] = imread(current_image);

    % Convert to Grayscale
    Gray_Map_A = rgb2gray(Image_A);
    %imshow(Gray_Map_A);
    Gray_Map_A = imgaussfilt(Gray_Map_A,5);
    %Click on Points
    figure(1)
    imshow(Gray_Map_A);
    % [x,y] = getpts;

    x = 1000 * [1.8808
                1.6585
                2.0770
                2.4213
                2.6001
                2.0813
                1.8242
                2.1947
                0.2332
                0.2637
                1.1006
                1.1049];

    y = 1000 * [0.8614
                1.0445
                1.7027
                1.6765
                2.3870
                2.2519
                2.7706
                2.9711
                3.4114
                3.0452
                1.4150
                2.4655];


    hold on
    plot(x,y,'*r')
    hold on

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %               Data Processing                 %
    %     (Finger Features Metric Processing)       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % ASSUMPTIONS:
    % 1- The image is divided in 6 features (5 fingers + knuckles)
    % 2- Each feature contains 2 points
    % 3- Points are selected in order of features 
    %       i.e., Points 1 and 2 correcpond to feature 1
    %             Points 3 and 4 correcpond to feature 2, etc..


    % Setting points 
    %   Col 1 = x axis
    %   Col 2 = y axis
    points = [x y];
    metrics = [];
    data = [];
    % Iterate through each feature
    for feature = 1:2:size(points,1)

        % TO-DO: Special scenario for knuckles
        %if feature > 10:
        
        % Calculations Strategy
        % 1- Find the line going through the 2 feature points
        % 2- Find the distance d between the 2 feature points
        % 3- Find the perpendicular at each point, and sweep
        %     the pixels stating at a distance of d/2 away
        % 4- Find the local minima and maxima of the swept points


        % Find Feature Vector going through both feature points
        pt_1 = points(feature,:);
        pt_2 = points(feature+1,:);

        %% Do the math
        % Get slope and y int of line AB
        slope = (pt_2(2)-pt_1(2)) / (pt_2(1)-pt_1(1)); 
        yint = pt_2(2) - slope*pt_2(1); 

        % Find distance between both feature points
        dist = pdist([pt_1; pt_2],'euclidean');

        % Slope of perpendiculare line
        perpSlope = -1 / slope;

        % For each point find the Feature Metric / Finger Thickness
        feature_metric = [];

        pts = [pt_1; pt_2];

        for i = 1:2
            mid_pt = pts(i,:);

            % Specify how wide of an area the sweep should do
            sweep_range = dist*0.67; 

            % If the feature is a finger
            if feature <= 10

                % Find the end points of the perpendicular line with length Clen*2
                sweep_dist = [(sweep_range*sqrt(1/(1+perpSlope^2))),... 
                              (perpSlope*sweep_range*sqrt(1/(1+perpSlope^2)))];
                    
                start_pt = mid_pt - sweep_dist;
                end_pt = mid_pt + sweep_dist;
            
            % If the feature is a knucle
            else

                % Find the end points of the perpendicular line with length Clen*2
                sweep_dist = [(dist*sqrt(1/(1+slope^2))),... 
                              (slope*dist*sqrt(1/(1+slope^2)))];


                mid_pt = (pt_1(:) + pt_2(:)).'/2
                start_pt = mid_pt - sweep_dist;
                end_pt = mid_pt + sweep_dist;

            end

            pt_plot = [start_pt; end_pt];
            pixel_data = improfile(Gray_Map_A,pt_plot(:,1),pt_plot(:,2),dist);

            curr_metric = 0;
            feature_metric = [feature_metric; curr_metric];

            figure(1)
            hold on
            plot(pt_plot(:,1),pt_plot(:,2),'--*m')
            hold on
        end 
        
        % Update metrics matrix
        metrics = [metrics, feature_metric];
    end
end    

t=1:1:size(pixel_data,1);
TF = islocalmax(pixel_data);
figure(2)
plot(pixel_data)
hold on
plot(t(TF),pixel_data(TF),'r*')
d=diff(pixel_data)
hold on