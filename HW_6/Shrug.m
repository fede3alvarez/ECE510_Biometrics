close all
clear all;
clc;

% Code will be available at
% https://github.com/fede3alvarez/ECE510_Biometrics
% Homework 6 - Shrug


%---------------------------------------
%  Video Selection and clean-up
%---------------------------------------

% Videos
 Available_videos = ["Video_00.MOV"];

% video_to_process = Available_videos(4,:)
% raw_video = VideoReader(video_to_process);
% firstFrame = read(raw_video,1);

% Note:
% Due to codecs issue in matlab, and matlab been a very sub-optimal tool,
% the developer had issue with the VideoReader tool (per internet forums
% the only fix is a OS update...). This issue is still present since 2017
% and even with matlab having a hefty price tag...
%
% To mitigate this issue, the developer manually decompose the video into
% frames, which will be used for this assigment.


%---------------------------------------
%  Read video frames
%---------------------------------------
Available_frames = ["Frames_Video_00/frame"];

temp_directory = "temp/";

% Choose 1st video
frames_to_process = Available_frames(1,:);

% Note:  
% In reality Video_00 = 1074 frames & Video_00 = 1054
num_of_frames = 1074;

%---------------------------------------
%  Global
%---------------------------------------

% Create a cascade detector object.
faceDetector = vision.CascadeObjectDetector();

% Should box dimensions as percentage of face dimensions
shoulder_height = 0.6;
shoulder_width = 0.5;
shoulder_start = 0.8;

% shrug Detection
shrug_detected = false;
shrug_message = "TBD";
shrug_frame_threshold = 5;
shrug_step = 5;
shrug_step_cnt = 0;
shrug_clear_cnt = 0;
shrug_clear = 45;
shoulders = zeros(num_of_frames,2);

%---------------------------------------
%  Iterate frame per frame
%---------------------------------------

for frame_num = 788:num_of_frames

    % Initialize frame specific variables
    %-------------------------------------
    rs_loc = [];   % Right shoulder points
    ls_loc = [];   % Left shoulder points


    % Load the current frame
    %------------------------------
    frame_path = strcat(frames_to_process,...
                        num2str(frame_num, '%04i'),...
                        '.png');

    [frame_Im, frame_Map] = imread(frame_path);
    frame_Im_gray = rgb2gray(frame_Im);
    frame_Im_raw = frame_Im;


%---------------------------------------
%  Locate head
%---------------------------------------

    % To locate a face, the following example is used as a model:
    % https://www.mathworks.com/help/vision/ug/face-detection-and-tracking-using-the-klt-algorithm.html

    % Detect the face
    face_box = step(faceDetector, frame_Im);

    % Decompose face_box for later
    % Assumption: 
    % IF multiple faced are found, the last one if the correct one
    face = size(face_box,1);
    face_start_x = face_box(face,1);
    face_start_y = face_box(face,2);
    face_width = face_box(face,3);
    face_height = face_box(face,4); 

    % Add face box to the image
    frame_Im = insertShape(frame_Im, "rectangle", face_box);

    %---------------------------------------
    %  Locate Shoulders Box
    %---------------------------------------

    % Right Shoulder
    %---------------------

    % Find area to sweep based on face_box
    rs_top_left = [face_start_x + face_width,...
                   face_start_y + round(shoulder_start * face_height)];

    rs_bottom_left  = [face_start_x + face_width,...
                      face_start_y +...
                           round((1+shoulder_height) * face_height)];

    rs_top_right = [face_start_x + round((1+shoulder_width) * face_width),...
                   face_start_y + round(shoulder_start * face_height)];

    rs_bottom_right= [face_start_x + round((1+shoulder_width) * face_width),...
                      face_start_y +...
                           round((1+shoulder_height) * face_height)];


    % Left Shoulder
    %---------------------


    % Find area to sweep based on face_box
    ls_top_left = [face_start_x - round(shoulder_width * face_width),...
                   face_start_y + round(shoulder_start * face_height)];

    ls_bottom_left = [face_start_x - round(shoulder_width * face_width),...
                      face_start_y +...
                           round((1+shoulder_height) * face_height)];

    ls_top_right = [face_start_x,...
                   face_start_y + round(shoulder_start * face_height)];

    ls_bottom_right = [face_start_x,...
                      face_start_y +...
                           round((1+shoulder_height) * face_height)];

    % Check we are not out of bounds...
    %-----------------------------------
    %  Out of bound through the bottom of the image
    if (rs_bottom_left(2) >= size(frame_Im_gray,1) ||...
        rs_bottom_right(2) >= size(frame_Im_gray,1) ||...
        ls_bottom_left(2) >= size(frame_Im_gray,1) ||...
        ls_bottom_right(2) >= size(frame_Im_gray,1))
       
            rs_bottom_left(2) = size(frame_Im_gray,1);
            rs_bottom_right(2) = size(frame_Im_gray,1);
            ls_bottom_left(2) = size(frame_Im_gray,1);
            ls_bottom_right(2) = size(frame_Im_gray,1);
    end

    %---------------------------------------
    %  Sweep should box &
    %     identify shoulders
    %---------------------------------------

    % Right Shoulder
    %---------------------

    % Iterate through "windows" until the box is swept
    for rs_pixel_x = rs_top_left(1):rs_top_right(1)

        % Get pixels intensity
        rs_int_y = frame_Im_gray(rs_top_left(2):rs_bottom_left(2),...
                                 rs_pixel_x);

        % Smooth out and get gradiant
        rs_int_y = smoothn(rs_int_y);
        rs_grad = gradient(rs_int_y);

        % Get min gradient point and record location
        [rs_min_val, rs_min_idx] = min(rs_grad);
        rs_loc = [rs_loc;...
                  [rs_pixel_x, (rs_top_left(2) + rs_min_idx)]];

    end % Sweep though all right shoulder slides


    coeff_rs_pts = polyfit(rs_loc(:,1),rs_loc(:,2),2);
    rs_pts_x = rs_loc(:,1);
    rs_pts_y = polyval(coeff_rs_pts, rs_pts_x);




    % Left Shoulder
    %---------------------

    % Iterate through "windows" until the box is swept
    for ls_pixel_x = ls_top_left(1):ls_top_right(1)

        % Get pixels intensity
        ls_int_y = frame_Im_gray(ls_top_left(2):ls_bottom_left(2),...
                                 ls_pixel_x);

        % Smooth out and get gradiant
        ls_int_y = smoothn(ls_int_y);
        ls_grad = gradient(ls_int_y);

        % Get min gradient point and record location
        [ls_min_val, ls_min_idx] = min(ls_grad);
        ls_loc = [ls_loc;...
                  [ls_pixel_x, (ls_top_left(2) + ls_min_idx)]];

    end % Sweep though all right shoulder slides


    coeff_ls_pts = polyfit(ls_loc(:,1),ls_loc(:,2),2);
    ls_pts_x = ls_loc(:,1);
    ls_pts_y = polyval(coeff_ls_pts, ls_pts_x);


    %---------------------------------------
    %  Process shoulders movement
    %    & identify shrug
    %---------------------------------------

    % Save shoulders position as a mean
    shoulders(frame_num,:) = [mean(rs_pts_y), mean(ls_pts_y)];
    
    % If the (mean for the last 5 frames) > (value 5 from frames ago)
    %   on both should shoulder, increase the counter
    if (frame_num > shrug_step)
        
        rs_curr_mean = mean(...
                      shoulders([(frame_num - shrug_step):frame_num],1));

        ls_curr_mean = mean(...
                      shoulders([(frame_num - shrug_step):frame_num],2));

        if(rs_curr_mean > shoulders((frame_num - shrug_step),1)) &...
          (ls_curr_mean > shoulders((frame_num - shrug_step),2))
            shrug_step_cnt = shrug_step_cnt + 1;
        else
            shrug_step_cnt = 0;
        end

    end

    % If the counter exceed the threshold, we have shrug!
    if (shrug_step_cnt >= shrug_frame_threshold)
        shrug_detected = true;
    end
    
    %---------------------------------------
    %  Plot
    %---------------------------------------

    % Data to plot Shoulder boxes
    rs_plot = [rs_top_left
               rs_top_right
               rs_bottom_right
               rs_bottom_left
               rs_top_left
              ];

    ls_plot = [ls_top_left
               ls_top_right
               ls_bottom_right
               ls_bottom_left
               ls_top_left
              ];

    % Plot Image
    f = figure('visible','off');
    %f = figure
    imshow(frame_Im)
    hold on
    plot(rs_plot(:,1), rs_plot(:,2), 'm-');
    hold on
    plot(ls_plot(:,1), ls_plot(:,2),'b-');
    hold on
    plot(rs_pts_x, rs_pts_y,'m*','LineWidth',2);
    hold on
    plot(ls_pts_x, ls_pts_y,'b*','LineWidth',2);


    processed_frame = strcat(temp_directory, ...
                             "Frame",...
                             num2str(frame_num, '%04i'),...
                             "_post_processing.png");

    if (shrug_detected)
        shrug_message = "SHRUG!";
        text(face_start_x + round(1.4*face_width),...
             face_start_y,...
             shrug_message,'Color','green','FontSize',60)

        % Keep the shrug detected message ON for 
        % shrug_clear frames
        shrug_clear_cnt = shrug_clear_cnt + 1;

        % If the messages was there too long
        %  clear everything and try again
        if(shrug_clear_cnt >= shrug_clear)
            shrug_detected = false;
            shrug_clear_cnt = 0;
            shrug_step_cnt = 0;
        end

    else
        shrug_message = "No shrug...";
         text(face_start_x + round(1.4*face_width),...
             face_start_y,...
             shrug_message,'Color','red','FontSize',60)
    end

    saveas(f,processed_frame)


end % End of frame per frame loop