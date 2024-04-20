close all
clear all;
clc;

% Code will be available at
% https://github.com/fede3alvarez/ECE510_Biometrics
% Homework 7 - Face Recognition


%---------------------------------------
%  Load Training Set
%---------------------------------------

% Load Training Set
dataset_folder = "faces_training";
dataset = imageDatastore(dataset_folder,...
                       "IncludeSubfolders", true,...
                       "LabelSource","foldernames");

% Split set between training and validation
%   Each folder / label / subject has ~20 pictures
[training_dataset, validation_dataset] = splitEachLabel(dataset, 7.0);

%[net,classNames] = imagePretrainedNetwork
net = googlenet
analyzeNetwork(net)