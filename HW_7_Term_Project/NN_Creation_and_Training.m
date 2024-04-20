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

% Split set between training (30%) and validation (70%)
%   Each folder / label / subject has ~20 pictures
[training_dataset, validation_dataset] = splitEachLabel(dataset, 7.0);


%---------------------------------------
%  Create Neural Network (GoogleNet)
%---------------------------------------

% Create Neural Network
%   The specific choice of Google Net is of not particular reason
%   Documentation and examples were able AND
%   Face Recognition is one of the targets of GoogleNet
neural_network = googlenet;
analyzeNetwork(neural_network);


% Update Neural Template for purpose (i.e. Face Recognition)

number_of_classes = numel(categories(training_dataset.Labels));

new_feature_learner = fullyConnectedLayer(number_of_classes, ...
                            "Name", "Facial Feature Learner",...
                            "WeightLearnRateFactor", 10,...
                            "BiasLearnRateFactor", 10);

new_classifier_layer = classificationLayer("Name", "Face Classifier");


% Update Layers
layer_graph = layerGraph(neural_network);

% Update Layer 142 to target Facial Features
feature_learner = neural_network.Layers(142);
layer_graph = replaceLayer(layer_graph, feature_learner.Name, new_feature_learner);

% Update Layer 144 to target Facial Classification
output_classifier = neural_network.Layers(144);
layer_graph = replaceLayer(layer_graph, output_classifier.Name, new_classifier_layer);

analyzeNetwork(layer_graph);


%---------------------------------------
%  Resize-up Image
%---------------------------------------
% Get Input size of 1st Layer
input_layer_size = neural_network.Layers(1).InputSize;

% Resize Image to fit within Neural Netwrok 
augmented_training_image = augmentedImageDatastore(...
                                        input_layer_size(1:2), ...
                                        training_dataset);

augmented_validation_image = augmentedImageDatastore(...
                                        input_layer_size(1:2), ...
                                        validation_dataset);


%---------------------------------------
%  Train Network 
%---------------------------------------

size_of_minibatch = 5;
validation_frequency = floor(...
           numel(augmented_training_image.Files) / size_of_minibatch...
                             );

training_options = trainingOptions("sgdm", ...
                        "MiniBatchSize", size_of_minibatch,...
                        "MaxEPochs", 6,...
                        "InitialLearnRate", 3e-4,...
                        "Shuffle", "every-epoch",...
                        "ValidationData", augmented_validation_image,...
                        "ValidationFrequency", validation_frequency,...
                        "Verbose", false,...
                        "Plots", "training-progress");

neural_network = trainNetwork(augmented_training_image, ...
                              layer_graph, ...
                              training_options);

% Save Neural Network for laater use
save neural_network;