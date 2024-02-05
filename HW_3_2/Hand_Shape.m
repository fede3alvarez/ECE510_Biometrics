close all
clear all;
clc;

% Code will be available at
% https://github.com/fede3alvarez/ECE510_Biometrics
% Homework 3 - Hand Shapes

% Images
Available_images = [ 'HandImage01.jpeg'
                    'HandImage02.jpeg'
                    'HandImage03.jpeg'
                    'HandImage04.jpeg'
                    'HandImage05.jpeg'
                    %'HandImage00.jpeg'
                    ];

Feature_Map = ['Pinky, Top    '
               'Pinky, Lower  '
               'Ring, Top     '
               'Ring, Lower   '
               'Middle, Top   '
               'Middle, Lower '
               'Index, Top    '
               'Index, Lower  '
               'Thumb, Top    '
               'Thumb, Lower  '
               'Knuckle, Right'
               'Knuckle, Left '
                ];

f_figure = 1;

for m = 1:size(Available_images,1)
    current_image = Available_images(m,:);
    [Image_A, Map_A] = imread(current_image);

    % Convert to Grayscale
    Gray_Map_A = rgb2gray(Image_A);
    
    % Filter Image
    %Gray_Map_A = imgaussfilt(Gray_Map_A,2);
    
    % Plot Image
    figure(f_figure)

    % % Fix Matlab axis so getpts is accurate
    % hax = axes('Parent', figure(f_figure));
    % axis(hax,'manual');
    imshow(Gray_Map_A);
    title(current_image);
        
    % % Get user to select Point
    % [x,y] = getpts;

    % % Save data
    % save_file_selections = strcat(current_image(1:11),'_data')
    % save(save_file_selections,'x','y')

    % Load saved data
    load_selections = strcat(current_image(1:11),'_data.mat');
    selections = importdata(load_selections);
    x = selections.x;
    y = selections.y;
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
    pixel_data = [];
    f_subplot = 1;

    % Iterate through each feature
    for feature = 1:2:size(points,1)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                 Data Processing               %
        %              (Image Processing to             %
        %             obtain feature metrics)           %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        pSlope = -1 / slope;

        % For each point find the Feature Metric / Finger Thickness
        feature_metric = [];

        pts = [pt_1; pt_2];

        for i = 1:2
            mid_pt = pts(i,:);

            % Specify how wide of an area the sweep should do
            % There is nothing special about this ratio - trial and error
            sweep_range = dist*0.67; 

            % If the feature is a finger
            if feature <= 10

                % Find the end points of the perpendicular line 
                %   with length sweep_range
                sweep_dist = [sweep_range*sqrt(1/(1+pSlope^2)),... 
                              pSlope*sweep_range*sqrt(1/(1+pSlope^2))];

                % Set Start & End points to sweep pixel data
                start_pt = mid_pt - sweep_dist;
                end_pt = mid_pt + sweep_dist;

            % If the feature is a knucle
            else

                % Find the end points of the perpendicular line with length Clen*2
                sweep_dist = [(dist*sqrt(1/(1+slope^2))),... 
                              (slope*dist*sqrt(1/(1+slope^2)))];

                % If the feature is a knuckle, 
                %   Plot / extend the sweep on the line marked by the 2
                %       feature points, and not on the perpendicular line
                %   Extend distance to sweep from the mid point (between
                %       the feature points) instead than from each point
                %       separetely (this will result in calculating the
                %       same distance / metric twice).
                mid_pt = (pt_1(:) + pt_2(:)).'/2;
                start_pt = mid_pt - sweep_dist;
                end_pt = mid_pt + sweep_dist;

            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 Data Processing               %
            %           At this point we have the           %
            %             obtain feature metrics)           %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Note:
            %   At this point we have the raw pixel data ready to be
            %       processed (i.e. pixel value of feature metric)
            %  This section involves obtaining the pixel data, and
            %       calculating its derivative. Based on the derivative,
            %       the feature metrics will be decided.

            pt_plot = [start_pt; end_pt];
            pixel_line = improfile(Gray_Map_A,pt_plot(:,1),pt_plot(:,2),dist);
            d_pixel = diff(pixel_line);
            [min_d_pixel_val, min_d_pixel_idx] = min(d_pixel); 
            [max_d_pixel_val, max_d_pixel_idx] = max(d_pixel); 

            curr_metric = abs(max_d_pixel_val - min_d_pixel_val);
            feature_metric = [feature_metric; curr_metric];

            % If the feature is a finger
            if feature <= 10

                % Find the end points of the perpendicular line 
                %   with length sweep_range
                metric_dist = [(curr_metric/2)*sqrt(1/(1+pSlope^2)),... 
                              pSlope*(curr_metric/2)*sqrt(1/(1+pSlope^2))];

                % Set Start & End points to sweep pixel data
                metric_start = mid_pt - metric_dist;
                metric_end = mid_pt + metric_dist;

            % If the feature is a knucle
            else

                % Find the end points of the perpendicular line with length Clen*2
                metric_dist = [((curr_metric/2)*sqrt(1/(1+slope^2))),... 
                              (slope*(curr_metric/2)*sqrt(1/(1+slope^2)))];

                % If the feature is a knuckle, 
                %   Plot / extend the sweep on the line marked by the 2
                %       feature points, and not on the perpendicular line
                %   Extend distance to sweep from the mid point (between
                %       the feature points) instead than from each point
                %       separetely (this will result in calculating the
                %       same distance / metric twice).
                mid_pt = (pt_1(:) + pt_2(:)).'/2;
                metric_start = mid_pt - metric_dist;
                metric_end = mid_pt + metric_dist;

            end

            metric_pts = [metric_start; metric_end];

            % Plot Hand Image with selected points
            %   AND sweep lines
            %   AND calculate feature metric
            figure(f_figure)
            hold on
            plot(pts(:,1),pts(:,2),'*b')
            hold on
            plot(pt_plot(:,1),pt_plot(:,2),'--*b')
            hold on
            plot(metric_pts(:,1),metric_pts(:,2),'*r:')
            hold on


            % For each feature in a separete figure, plot
            %   The pixel values for the line to be computed
            %   The pixel dereivatives, and calculated max / mins
            %       (max/mins correspond to beginning and end of metrics)
            figure((f_figure+1))
            hold on
            subplot(6,2,f_subplot)
            plot(pixel_line,'b*:')
            hold on
            plot(d_pixel,'g+:')
            hold on
            plot(min_d_pixel_idx,d_pixel(min_d_pixel_idx),'rd')
            hold on
            plot(max_d_pixel_idx,d_pixel(max_d_pixel_idx),'rd')
            grid on
            subplot_title = strcat(current_image(10:11),', ',...
                            Feature_Map(feature,:))
            title(subplot_title);
            f_subplot = f_subplot+1;


        end     % All points in a feature - Iteration Done

        % Update metrics matrix
        metrics = [metrics, feature_metric];
    
        
    hold on
    
    end     % All features Iteration
    f_figure = f_figure + 2;
end         % All images Iteration

% t=1:1:size(pixel_data,1);
% TF = islocalmax(pixel_data);
% figure(2)
% plot(pixel_data)
% hold on
% plot(t(TF),pixel_data(TF),'r*')
% d=diff(pixel_data)
% hold on