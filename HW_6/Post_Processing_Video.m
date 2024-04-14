close all
clear all;
clc;

% Code will be available at
% https://github.com/fede3alvarez/ECE510_Biometrics
% Homework 6 - Shrug

% This video was written using the following guide:
% https://www.mathworks.com/help/matlab/import_export/convert-between-image-sequences-and-video.html

% Directory where to collect the frames
temp_directory = "temp/";
frame_extension = ".png";
frame_rate = 60;

% Matlab magic to read all frames
%processed_frames = dir(fullfile(temp_directory,"Frames*"));
frame_search = strcat(temp_directory,"Frame*",frame_extension);
processed_frames = dir(frame_search)

% Video
outputVideo = VideoWriter("Video_00_Shrug.avi");
outputVideo.FrameRate = frame_rate;
open(outputVideo)


for frame = 1:length(processed_frames)
   img = imread(strcat(temp_directory,processed_frames(frame).name));
   writeVideo(outputVideo,img)
end


close(outputVideo)