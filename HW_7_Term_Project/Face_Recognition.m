close all
clear all;
clc;

% Load Neural Network
load neural_network;

% Folder where Test Images are located
test_images_folder = 'Test_Images';

% Iterate through all test images
for test_image = 1:3

    % Iterate through all test images
    image = strcat(test_images_folder,...
                   "/Test_",...
                   num2str(test_image), ...
                   '.jpg');

    % Test image through neural networks
    test_network(neural_network, image);
end


% Process Image Through Neural Network
function test_network(net, image)

    % Read and resize image to fit with
    %   The NN expected input size
    I = imread(image);
    G = imresize(I, [224, 224]);
    
    % Classify Test Images
    [Label, Prob] = classify(net, G);

    % Plot Image with Classification and Prob as Title
    figure;
    imshow(G);
    title(strcat(char(Label),'-', ...
                 num2str(max(Prob)*100,'%4.2f'), '%'));
    pause(1);
end