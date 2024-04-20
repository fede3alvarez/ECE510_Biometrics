close all
clear all;
clc;

% Code will be available at
% https://github.com/fede3alvarez/ECE510_Biometrics
% Homework 7 - Face Recognition


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
image = "Sample_Face.png";


% Note:  
% In reality Video_00 = 1074 frames & Video_00 = 1054
num_of_frames = 1074;

%---------------------------------------
%  Global
%---------------------------------------

% Create a cascade detector object.
faceDetector = vision.CascadeObjectDetector();


%---------------------------------------
%  Iterate frame per frame
%---------------------------------------


[frame_Im, frame_Map] = imread(image);
frame_Im_gray = rgb2gray(frame_Im);
frame_Im_raw = frame_Im;


%---------------------------------------
%  Locate head
%---------------------------------------

% To locate a face, the following example is used as a model:
% https://www.mathworks.com/help/vision/ug/face-detection-and-tracking-using-the-klt-algorithm.html

% Detect the face
face_box = step(faceDetector, frame_Im);

% Add face box to the image
frame_Im = insertShape(frame_Im, "rectangle", face_box);


    
%---------------------------------------
%  Plot
%---------------------------------------

% Plot Image
%f = figure('visible','off');
f = figure
imshow(frame_Im)
hold on
