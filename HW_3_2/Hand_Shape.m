close all
clear all;
clc;

% Code will be available at
% https://github.com/fede3alvarez/ECE510_Biometrics
% Homework 3 - Hand Shapes

% Images
Available_images = ['HandImage01.jpeg'
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

% For each point find the Feature Metric / Finger Thickness
feature_metric = [];

for m = 1:size(Available_images,1)
    current_image = Available_images(m,:);
    [Image_A, Map_A] = imread(current_image);

    % Convert to Grayscale
    Gray_Map_A = rgb2gray(Image_A);
    
    % Filter Image
    %Gray_Map_A = imgaussfilt(Gray_Map_A,6);
    
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
        user_selected_pt_1 = points(feature,:);
        user_selected_pt_2 = points(feature+1,:);

        %% Do the math
        % Get slope and y int of line AB
        slope = (user_selected_pt_2(2)-user_selected_pt_1(2)) / ...
                (user_selected_pt_2(1)-user_selected_pt_1(1)); 
        yint = user_selected_pt_2(2) - slope*user_selected_pt_2(1); 

        % Find distance between both feature points
        dist = pdist([user_selected_pt_1;user_selected_pt_2],'euclidean');

        % Slope of perpendiculare line
        pSlope = -1 / slope;

        user_selected_pts = [user_selected_pt_1; user_selected_pt_2];

        for i = 1:2
            mid_pt = user_selected_pts(i,:);

            % Specify how wide of an area the sweep should do
            % There is nothing special about this ratio - trial and error
            sweep_range = dist*0.67; 

            % If the feature is a finger
            if feature <= 10

                % Specify how wide of an area the sweep should do
                % There is nothing special about this ratio - trial and error
                sweep_range = dist*0.67; 

                % Find the end points of the perpendicular line 
                %   with length sweep_range
                sweep_dist = [sweep_range*sqrt(1/(1+pSlope^2)),... 
                              pSlope*sweep_range*sqrt(1/(1+pSlope^2))];

                % Set Start & End points to sweep pixel data
                start_pt = mid_pt - sweep_dist;
                end_pt = mid_pt + sweep_dist;

            % If the feature is a knucle
            else

                % Specify how wide of an area the sweep should do
                % There is nothing special about this ratio - trial and error
                sweep_range = dist*1.4; 

                % Find the end points of the perpendicular line 
                sweep_dist = [(sweep_range*sqrt(1/(1+slope^2))),... 
                              (slope*sweep_range*sqrt(1/(1+slope^2)))];

                % If the feature is a knuckle, 
                %   Plot / extend the sweep on the line marked by the 2
                %       feature points, and not on the perpendicular line
                %   Extend distance to sweep from the mid point (between
                %       the feature points) instead than from each point
                %       separetely (this will result in calculating the
                %       same distance / metric twice).
                mid_pt = (user_selected_pt_1(:) +...
                          user_selected_pt_2(:)).'/2;
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

            sweep_pts = [start_pt; end_pt];
            pixel_line = improfile(Gray_Map_A,...
                                   sweep_pts(:,1),...
                                   sweep_pts(:,2));
            d_pixel = diff(pixel_line);
            [min_val, min_idx] = min(d_pixel); 
            [max_val, max_idx] = max(d_pixel);
            
            d2_pixel = diff(d_pixel);
            [min2_val, min2_idx] = min(d2_pixel); 
            [max2_val, max2_idx] = max(d2_pixel);

            curr_metric = pdist([max_idx; min_idx],'euclidean');
            feature_metric = [feature_metric; curr_metric];

            f_plot = [min_idx, pixel_line(min_idx)
                      max_idx, pixel_line(min_idx)];

            figure(f_figure)
            hold on
            plot(user_selected_pts(:,1),user_selected_pts(:,2),'*b')
            hold on

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

                % Find the end points of the perpendicular line 
                metric_dist = [((curr_metric/2)*sqrt(1/(1+slope^2))),... 
                              (slope*(curr_metric/2)*sqrt(1/(1+slope^2)))];

                % If the feature is a knuckle, 
                %   Plot / extend the sweep on the line marked by the 2
                %       feature points, and not on the perpendicular line
                %   Extend distance to sweep from the mid point (between
                %       the feature points) instead than from each point
                %       separetely (this will result in calculating the
                %       same distance / metric twice).
                mid_pt = (user_selected_pt_1(:) + ...
                          user_selected_pt_2(:)).'/2;
                user_selected_pt_1 = mid_pt;
                user_selected_pt_2 = mid_pt;
                metric_start = mid_pt - metric_dist;
                metric_end = mid_pt + metric_dist;

            end

            metric_pts = [metric_start; metric_end];

            % Plot Hand Image with selected points
            %   AND sweep lines
            %   AND calculate feature metric
            figure(f_figure)
            hold on
            plot(user_selected_pts(:,1),user_selected_pts(:,2),'*w')
            hold on
            plot(sweep_pts(:,1),sweep_pts(:,2),'--*b')
            hold on
            plot(metric_start(:,1),metric_start(:,2),'dr:')
            hold on
            plot(metric_end(:,1),metric_end(:,2),'r:square')
            hold on
            plot(metric_pts(:,1),metric_pts(:,2),'r:')
            hold on
            plot(mid_pt(:,1),mid_pt(:,2),'g:o')
            hold on
            legend('Points Selected',...
                   'Area to be sweep / analyzed',...
                   'Calculated Feature Metric')


            % For each feature in a separete figure, plot
            %   The pixel values for the line to be computed
            %   The pixel dereivatives, and calculated max / mins
            %       (max/mins correspond to beginning and end of metrics)
            figure((f_figure+1))
            hold on
            subplot(6,2,f_subplot)
            yyaxis left
            plot(pixel_line,'b*:')
            hold on
            plot(f_plot(:,1),f_plot(:,2),'msquare-','LineWidth',3)
            hold on
            yyaxis right
            plot(d_pixel,'g+:')
            hold on
            plot(min_idx,d_pixel(min_idx),'rd','LineWidth',3)
            hold on
            plot(max_idx,d_pixel(max_idx),'rd','LineWidth',3)
            hold on
            % plot(d2_pixel,'msquare:')
            % hold on
            % plot(min2_idx,d2_pixel(min2_idx),'ko')
            % hold on
            % plot(max2_idx,d2_pixel(max2_idx),'ko')
            grid on
            legend('Area Scanned',...
                   'Calculated Feature Lenght',...
                   'Derivative',...
                   'Min/Max Der Value')
            subplot_title = strcat(current_image(10:11),', ',...
                            Feature_Map(f_subplot,:));
            title(subplot_title);
            f_subplot = f_subplot+1;


        end     % All points in a feature - Iteration Done

    
        
    hold on
    
    end     % All features Iteration
    f_figure = f_figure + 2;
end         % All images Iteration


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Pairwise comparison              %
%            among the 5 hand images            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Re-arrange the metrics collects on a per image (cols) basis
feature_matrix = reshape(feature_metric,12,5);

% Create a Verification / Pairwase Comparison Matrix
verif_matrix = zeros(5,5);

% Fill each entry with the corresponding comparison
for feature_to_compare = 1:size(Available_images,1)
    for comparable = (feature_to_compare+1):size(Available_images,1)
        
        % Perform the pairwise euclidean comparison
        verif_matrix(feature_to_compare,comparable) = ...
             pdist(...
                   [feature_matrix(:,feature_to_compare)';...
                    feature_matrix(:,comparable)'],...
                   'euclidean');

    end % No more comaprables
end % No more features to compare

% Display the Matrix
verif_matrix