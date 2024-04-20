close all
clear all;
clc;


load neural_network;

test_images_folder = 'Test_Images';

for test_image = 1:3
    image = strcat(test_images_folder,...
                   "/Test_",...
                   num2str(test_image), ...
                   '.jpg');
    test_network(neural_network, image)
end


function test_network(net, image)

    I = imread(image);
    G = imresize(I, [224, 224]);
    
    [Label, Prob] = classify(net, G)
    max(Prob)
    num2str(max(Prob),2)

    f = figure;
    imshow(G);
    title(strcat(char(Label),'-', ...
                 num2str(max(Prob)*100,'%4.2f'), '%'));
    pause(1);
end