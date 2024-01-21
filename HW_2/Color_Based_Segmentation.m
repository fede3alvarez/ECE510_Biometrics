close all
clear all;
clc;

% Code will be available at
% https://github.com/fede3alvarez/ECE510_Biometrics

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Parameters of interest                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k_values = [2 3 5]; % Array with the numbers of k-clusters to test
k_repeats = 100;    % # of repeat of the kmeans for each k cluster value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Read Image (*.jpg) and massage it       %
%             to make matlab happy...            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%/home/fico/Repos/ECE510_Biometrics/HW_2/Images/00B_Hand.jpg\

% Note: jpeg images are read directly as RGB  matrices
%    https://www.mathworks.com/matlabcentral/answers/23119-jpg-to-rgb-image
[Image_A, Map_A] = imread('00B_Hand.jpg');
Image = double(Image_A);

% Separete data in colors
[rows, cols, colors]  = size(Image);
red = Image_A(:, :, 1);
green = Image_A(:, :, 2);
blue = Image_A(:, :, 3);

% Rearrange data for kmeans clustering
data = double([red(:), green(:), blue(:)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              K-Means Clustering                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In this section we setup the numbers of K cluster to be used

% TO-DO: Implement this as a for-loop (2,3,5)
k = 2;
% Prob_Map_k_2
% Plot_Map_2

% From: Homework Prompt:
% For example, for k=3, and 100 repetitions, the output should be 3 maps:
%       P(x belongs to cluster i) | {Rx,Gx,Bx} ), for i={1,2,3} 
%       Note: ({Rx,Gx,Bx} are the chromatic values of x). 
%
% Desing: This will be implement as k-probability maps initialized to zero
%         The k-prob matrix will have dimensions rows x cols x k
%         Each time a pixel is assigned to a cluster, the corresponding
%         map entry will be increased by 1/k_repeats
Prob_Map = zeros(rows,cols,k);

%Repeat this for 100
for i = 1:k_repeats 

    % Setting random k initial conditions
    init_cond = [];
    for j = 1:k 
        r = randi(size(data,1),1);
        init_cond = [init_cond; data(r,:)];
    end

    % Sort clusters based on Red, and the Green values (Higher 1st / Top)
    init_cond = sortrows(init_cond,'descend');

    % Run k-means
    [idx, C] = kmeans(data, k,'start', init_cond);

    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %              Probability Matrix                %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %
    % Populate Probability Matrix
    % The k-prob matrix has dimensions rows x cols x k
    % We populate it when iterating through each k clusters
    % and adding on each iteration a value of 1/(k_repeats * k)
    % to the assigned cluster
    %
    %   Notes: 
    %    * 1 = k_repeats (total testing times) * 1/k_repeats (each test)
    %    * The 1/k factor is aded to correct the fact that Matlab returns
    %        the index of cluster and not a 0 OR 1 value.
    %    * Clusters were all ready sorted at the inital conditions, before
    %        running kmeans,so Matlab is returning an already sorted result
    for j = 1:k 
        Prob_Map(:,:,j) = Prob_Map(:,:,j)...
                        +((1/(j*k_repeats)) * reshape(idx == j,rows,cols));
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Plots                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

% Create a plot map based on the Probability Map
Plot_Map = ones(rows,cols);

% Plot_Map starts as a matrix of zeros (for cluster 1 as default)
% For each subsequent cluster, 
%   if (prob map value of cluster j) > (prob map value of cluster j-1)
%   then, update the entry of plot map to the corresponding cluster
Plot_Map(Prob_Map(:,:,j) > Plot_Map) = j;
for j = 2:k 
        Plot_Map(Prob_Map(:,:,(j)) > Prob_Map(:,:,(j-1))) = j;
end


% Analysis of inputs for debugging
figure(1)
subplot(2, 2, 1);
imshow(Image_A);
title('Original Color Image');
subplot(2, 2, 2);
imshow(red);
title('Red Image');
subplot(2, 2, 3);
imshow(green);
title('Green Image');
subplot(2, 2, 4);
imshow(blue);
title('Blue Image');

% Plotting of clusters and results
figure(2)
imshow(Plot_Map);