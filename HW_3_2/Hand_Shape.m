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
                    'HandImage00.jpeg'
                    ];

f_figure = 1;

for m = 1:size(Available_images,1)
    current_image = Available_images(m,:);
    [Image_A, Map_A] = imread(current_image);

    % Convert to Grayscale
    Gray_Map_A = rgb2gray(Image_A);
    
    % Filter Image
    %Gray_Map_A = imgaussfilt(Gray_Map_A,5);
    
    % Plot Image
    figure(f_figure)

    % Fix Matlab axis so getpts is accurate
    % hax = axes('Parent', figure(f_figure));
    % axis(hax,'manual');
    imshow(Gray_Map_A);
        
    % Get user to select Point
    % [x,y] = getpts;
    % 
    % % Save data
    % save_file_selections = strcat(current_image(1:11),'_data')
    % save(save_file_selections,'x','y')

    % Load saved data
    load_selections = strcat(current_image(1:11),'_data.mat')
    selections = importdata(load_selections)
    hold on

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %               Data Processing                 %
    % %     (Finger Features Metric Processing)       %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    % 
    % % ASSUMPTIONS:
    % % 1- The image is divided in 6 features (5 fingers + knuckles)
    % % 2- Each feature contains 2 points
    % % 3- Points are selected in order of features 
    % %       i.e., Points 1 and 2 correcpond to feature 1
    % %             Points 3 and 4 correcpond to feature 2, etc..
    % 
    % 
    % % Setting points 
    % %   Col 1 = x axis
    % %   Col 2 = y axis
    % points = [selections.x selections.y];
    % metrics = [];
    % data = [];
    % pixel_data = [];
    % f_subplot = 1;
    % 
    % % Iterate through each feature
    % for feature = 1:2:size(points,1)
    % 
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     %                 Data Processing               %
    %     %              (Image Processing to             %
    %     %             obtain feature metrics)           %
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    %     % Calculations Strategy
    %     % 1- Find the line going through the 2 feature points
    %     % 2- Find the distance d between the 2 feature points
    %     % 3- Find the perpendicular at each point, and sweep
    %     %     the pixels stating at a distance of d/2 away
    %     % 4- Find the local minima and maxima of the swept points
    % 
    % 
    %     % Find Feature Vector going through both feature points
    %     pt_1 = points(feature,:);
    %     pt_2 = points(feature+1,:);
    % 
    %     %% Do the math
    %     % Get slope and y int of line AB
    %     slope = (pt_2(2)-pt_1(2)) / (pt_2(1)-pt_1(1)); 
    %     yint = pt_2(2) - slope*pt_2(1); 
    % 
    %     % Find distance between both feature points
    %     dist = pdist([pt_1; pt_2],'euclidean');
    % 
    %     % Slope of perpendiculare line
    %     pSlope = -1 / slope;
    % 
    %     % For each point find the Feature Metric / Finger Thickness
    %     feature_metric = [];
    % 
    %     pts = [pt_1; pt_2];
    % 
    %     for i = 1:2
    %         mid_pt = pts(i,:);
    % 
    %         % Specify how wide of an area the sweep should do
    %         % There is nothing special about this ratio - trial and error
    %         sweep_range = dist*0.67; 
    % 
    %         % If the feature is a finger
    %         if feature <= 10
    % 
    %             % Find the end points of the perpendicular line 
    %             %   with length sweep_range
    %             sweep_dist = [sweep_range*sqrt(1/(1+pSlope^2)),... 
    %                           pSlope*sweep_range*sqrt(1/(1+pSlope^2))];
    % 
    %             % Set Start & End points to sweep pixel data
    %             start_pt = mid_pt - sweep_dist;
    %             end_pt = mid_pt + sweep_dist;
    % 
    %         % If the feature is a knucle
    %         else
    % 
    %             % Find the end points of the perpendicular line with length Clen*2
    %             sweep_dist = [(dist*sqrt(1/(1+slope^2))),... 
    %                           (slope*dist*sqrt(1/(1+slope^2)))];
    % 
    %             % If the feature is a knuckle, 
    %             %   Plot / extend the sweep on the line marked by the 2
    %             %       feature points, and not on the perpendicular line
    %             %   Extend distance to sweep from the mid point (between
    %             %       the feature points) instead than from each point
    %             %       separetely (this will result in calculating the
    %             %       same distance / metric twice).
    %             mid_pt = (pt_1(:) + pt_2(:)).'/2;
    %             start_pt = mid_pt - sweep_dist;
    %             end_pt = mid_pt + sweep_dist;
    % 
    %         end
    % 
    %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         %                 Data Processing               %
    %         %           At this point we have the            %
    %         %             obtain feature metrics)           %
    %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    %         % Note:
    %         %   At this point we have the raw pixel data ready to be
    %         %       processed (i.e. pixel value of feature metric)
    %         %  This section involves obtaining the pixel data, and
    %         %       calculating its derivative. Based on the derivative,
    %         %       the feature metrics will be decided.
    % 
    %         pt_plot = [start_pt; end_pt];
    %         pixel_line = improfile(Gray_Map_A,pt_plot(:,1),pt_plot(:,2),dist);
    %         d_pixel = diff(pixel_line);
    % 
    %         curr_metric = 0;
    %         feature_metric = [feature_metric; curr_metric];
    % 
    %         figure(f_figure)
    %         hold on
    %         plot(pt_plot(:,1),pt_plot(:,2),'--*m')
    %         hold on
    % 
    %         figure((f_figure+1))
    %         hold on
    %         subplot(6,2,f_subplot)
    %         plot(pixel_line,'b*:')
    %         hold on
    %         plot(d_pixel,'g+:')
    %         grid on
    %         f_subplot = f_subplot+1;
    % 
    % 
    %     end % All points in a feature Iteration
    % 
    %     % Update metrics matrix
    %     metrics = [metrics, feature_metric];
    % end    
    % plot(selections.x,selections.y,'*r')
    % hold on
    % 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %               Data Processing                 %
    % %     (Finger Features Metric Processing)       %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    % 
    % % ASSUMPTIONS:
    % % 1- The image is divided in 6 features (5 fingers + knuckles)
    % % 2- Each feature contains 2 points
    % % 3- Points are selected in order of features 
    % %       i.e., Points 1 and 2 correcpond to feature 1
    % %             Points 3 and 4 correcpond to feature 2, etc..
    % 
    % 
    % % Setting points 
    % %   Col 1 = x axis
    % %   Col 2 = y axis
    % points = [selections.x selections.y];
    % metrics = [];
    % data = [];
    % pixel_data = [];
    % f_subplot = 1;
    % 
    % % Iterate through each feature
    % for feature = 1:2:size(points,1)
    % 
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     %                 Data Processing               %
    %     %              (Image Processing to             %
    %     %             obtain feature metrics)           %
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    %     % Calculations Strategy
    %     % 1- Find the line going through the 2 feature points
    %     % 2- Find the distance d between the 2 feature points
    %     % 3- Find the perpendicular at each point, and sweep
    %     %     the pixels stating at a distance of d/2 away
    %     % 4- Find the local minima and maxima of the swept points
    % 
    % 
    %     % Find Feature Vector going through both feature points
    %     pt_1 = points(feature,:);
    %     pt_2 = points(feature+1,:);
    % 
    %     %% Do the math
    %     % Get slope and y int of line AB
    %     slope = (pt_2(2)-pt_1(2)) / (pt_2(1)-pt_1(1)); 
    %     yint = pt_2(2) - slope*pt_2(1); 
    % 
    %     % Find distance between both feature points
    %     dist = pdist([pt_1; pt_2],'euclidean');
    % 
    %     % Slope of perpendiculare line
    %     pSlope = -1 / slope;
    % 
    %     % For each point find the Feature Metric / Finger Thickness
    %     feature_metric = [];
    % 
    %     pts = [pt_1; pt_2];
    % 
    %     for i = 1:2
    %         mid_pt = pts(i,:);
    % 
    %         % Specify how wide of an area the sweep should do
    %         % There is nothing special about this ratio - trial and error
    %         sweep_range = dist*0.67; 
    % 
    %         % If the feature is a finger
    %         if feature <= 10
    % 
    %             % Find the end points of the perpendicular line 
    %             %   with length sweep_range
    %             sweep_dist = [sweep_range*sqrt(1/(1+pSlope^2)),... 
    %                           pSlope*sweep_range*sqrt(1/(1+pSlope^2))];
    % 
    %             % Set Start & End points to sweep pixel data
    %             start_pt = mid_pt - sweep_dist;
    %             end_pt = mid_pt + sweep_dist;
    % 
    %         % If the feature is a knucle
    %         else
    % 
    %             % Find the end points of the perpendicular line with length Clen*2
    %             sweep_dist = [(dist*sqrt(1/(1+slope^2))),... 
    %                           (slope*dist*sqrt(1/(1+slope^2)))];
    % 
    %             % If the feature is a knuckle, 
    %             %   Plot / extend the sweep on the line marked by the 2
    %             %       feature points, and not on the perpendicular line
    %             %   Extend distance to sweep from the mid point (between
    %             %       the feature points) instead than from each point
    %             %       separetely (this will result in calculating the
    %             %       same distance / metric twice).
    %             mid_pt = (pt_1(:) + pt_2(:)).'/2;
    %             start_pt = mid_pt - sweep_dist;
    %             end_pt = mid_pt + sweep_dist;
    % 
    %         end
    % 
    %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         %                 Data Processing               %
    %         %           At this point we have the            %
    %         %             obtain feature metrics)           %
    %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    %         % Note:
    %         %   At this point we have the raw pixel data ready to be
    %         %       processed (i.e. pixel value of feature metric)
    %         %  This section involves obtaining the pixel data, and
    %         %       calculating its derivative. Based on the derivative,
    %         %       the feature metrics will be decided.
    % 
    %         pt_plot = [start_pt; end_pt];
    %         pixel_line = improfile(Gray_Map_A,pt_plot(:,1),pt_plot(:,2),dist);
    %         d_pixel = diff(pixel_line);
    % 
    %         curr_metric = 0;
    %         feature_metric = [feature_metric; curr_metric];
    % 
    %         figure(f_figure)
    %         hold on
    %         plot(pt_plot(:,1),pt_plot(:,2),'--*m')
    %         hold on
    % 
    %         figure((f_figure+1))
    %         hold on
    %         subplot(6,2,f_subplot)
    %         plot(pixel_line,'b*:')
    %         hold on
    %         plot(d_pixel,'g+:')
    %         grid on
    %         f_subplot = f_subplot+1;
    % 
    % 
    %     end % All points in a feature Iteration
    % 
    %     % Update metrics matrix
    %     metrics = [metrics, feature_metric];
    % end     % All features Iteration
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