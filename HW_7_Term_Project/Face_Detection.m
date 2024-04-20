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


%---------------------------------------
%  Create Neural Network (GoogleNet)
%---------------------------------------

% Create Neural Network
net = googlenet;
analyzeNetwork(net)

% Get Input size of 1st Layer
Input_Layer_Size = net.Layers(1).InputSize;

% Update Layers as needed
Layer_Graph = layerGraph(net);

Feature_Learner = net.Layers(142);
Output_Classifier = net.Layers(144);

Number_of_Classes = numel(categories(training_dataset.Labels));

New_Feature_Learner = fullyConnectedLayer(Number_of_Classes, ...
                            "Name", "Facial Feature Learner",...
                            "WeightLearnRateFactor", 10,...
                            "BiasLearnRateFactor", 10);

New_Classifier_Layer = classificationLayer("Name", "Face Classifier");

Layer_Graph = replaceLayer(Layer_Graph, Feature_Learner.Name, New_Feature_Learner);
Layer_Graph = replaceLayer(Layer_Graph, Output_Classifier.Name, New_Classifier_Layer);

analyzeNetwork(Layer_Graph);

%---------------------------------------
%  Clean-up Image
%---------------------------------------

Pixel_Range = [-30 30];
Scale_Range = [0.9 1.1];

Image_Augmented = imageDataAugmenter('RandXReflection', true,...
                                   'RandXTranslation', Pixel_Range,...
                                   'RandYTranslation', Pixel_Range,...
                                   'RandXScale', Scale_Range,...
                                   'RandyScale', Scale_Range);

% Resize Image
Augmented_Training_Image = augmentedImageDatastore(...
                                        Input_Layer_Size(1:2), ...
                                        training_dataset,...
                                        'DataAugmentation',...
                                        Image_Augmented);

Augmented_Validation_Image = augmentedImageDatastore(...
                                        Input_Layer_Size(1:2), ...
                                        validation_dataset);

%---------------------------------------
%  Train Network 
%---------------------------------------

Size_of_Minibatch = 5;
Validation_Frequency = floor(numel(Augmented_Training_Image.Files) / Size_of_Minibatch);

Training_Options = trainingOptions("sgdm", ...
                        "MiniBatchSize", Size_of_Minibatch,...
                        "MaxEPochs", 6,...
                        "InitialLearnRate", 3e-4,...
                        "Shuffle", "every-epoch",...
                        "ValidationData", Augmented_Validation_Image,...
                        "ValidationFrequency", Validation_Frequency,...
                        "Verbose", false,...
                        "Plots", "training-progress");

net = trainNetwork(Augmented_Training_Image, Layer_Graph, Training_Options);