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

color_map = [255 0   0    % Red
             0   255 0    % Blue
             0   0   255  % Green
             0   0   0    % Black
             255 255 255  % White
             255 0   255  % Pink
             255 255 0];  % Brown

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


for k = [2, 3 5]
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
    
    
    % Repeat this loop for 100
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
    % Find the maximun probability for each cluster
    max_prob = max(Prob_Map,[],3);
    
    red_cluster_results = zeros(rows,cols);
    blue_cluster_results = zeros(rows,cols);
    green_cluster_results = zeros(rows,cols);
    
    for j = 1:k
        idx = (Prob_Map(:,:,j) == max_prob);
        
        red_cluster_results = red_cluster_results + idx*color_map(j,1);
        green_cluster_results = green_cluster_results + idx*color_map(j,2);
        blue_cluster_results = blue_cluster_results + idx*color_map(j,3);
    end
    
    Plot_Map = cat(3, ...
               red_cluster_results,...
               green_cluster_results,...
               blue_cluster_results);

    % Save data - need to find a more elegant way to do this
    if k == 2
        Prob_Map_k_2 = Prob_Map;
        Plot_Map_2 = Plot_Map;
    elseif k == 3
        Prob_Map_k_3 = Prob_Map;
        Plot_Map_3 = Plot_Map;
    elseif k == 5
        Prob_Map_k_5 = Prob_Map;
        Plot_Map_5 = Plot_Map;
    end
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
image(Plot_Map_2);
title('2 Clusters');

figure(3)
image(Plot_Map_3);
title('3 Clusters');

figure(4)
image(Plot_Map_5);
title('5 Clusters');
