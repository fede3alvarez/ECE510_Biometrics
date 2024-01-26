close all
clear all;
clc;

% Code will be available at
% https://github.com/fede3alvarez/ECE510_Biometrics

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Parameters of interest                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k_values = [2 3 5];       % Array with the # of k-clusters to test
k_repeats = 100;          % # of times to repeat kmeans
prob_limit = 0;

color_map = [255 0   0    % Red
             0   255 0    % Blue
             0   0   255  % Green
             0   0   0    % Black
             255 0   255  % Pink
             255 255 0];  % Brown

p = 1;                    % To keep track / initialize of plots / figures

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Read Image (*.jpg) and massage it       %
%             to make matlab happy...            %
%                 also plot it! 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Images
Available_images = [%'Hand00.jpg'
                    'Hand01.jpg'
                    %'Hand02.jpg'
                    'Hand03.jpg'
                    ];

for m = 1:size(Available_images,1)
    current_image = Available_images(m,:);
    
    % Note: jpeg images are read directly as RGB  matrices
    %    https://www.mathworks.com/matlabcentral/answers/23119-jpg-to-rgb-image
    [Image_A, Map_A] = imread(current_image);
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
    
    for k = [2, 3, 5]
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

        % Remove any probability below ther probability limnit
        Prob_Map(Prob_Map<prob_limit) = 0;            

        
        % Save data - need to find a more elegant way to do this
        if k == 2
            Prob_Map_k_2 = Prob_Map;
        
            % Plotting of Probability Map
            figure(p)
            h = heatmap(Prob_Map(:,:,1),'ColorLimits',[0 1]);
            h.Title = strcat('k=2 - Cluster 1 Prob - ',current_image);
            p = p + 1;
            
            figure(p)
            h = heatmap(Prob_Map(:,:,2),'ColorLimits',[0 1]);
            h.Title = strcat('k=2 - Cluster 2 Prob - ',current_image);
            p = p + 1;
        
        elseif k == 3
            Prob_Map_k_3 = Prob_Map;
        
            % Plotting of Probability Map
            figure(p)
            h = heatmap(Prob_Map(:,:,1),'ColorLimits',[0 1]);
            h.Title = strcat('k=3 - Cluster 1 Prob - ',current_image);
            p = p + 1;
            
            figure(p)
            h = heatmap(Prob_Map(:,:,2),'ColorLimits',[0 1]);
            h.Title = strcat('k=3 - Cluster 2 Prob - ',current_image);
            p = p + 1;
            
            figure(p)
            h = heatmap(Prob_Map(:,:,3),'ColorLimits',[0 1]);
            h.Title = strcat('k=3 - Cluster 3 Prob - ',current_image);
            p = p + 1;

        elseif k == 5
            
            Prob_Map_k_5 = Prob_Map;
       
            % Plotting of Probability Map
            figure(p)
            h = heatmap(Prob_Map(:,:,1),'ColorLimits',[0 1]);
            h.Title = strcat('k=5 - Cluster 1 Prob - ',current_image);
            p = p + 1;
            
            figure(p)
            h = heatmap(Prob_Map(:,:,2),'ColorLimits',[0 1]);
            h.Title = strcat('k=5 - Cluster 2 Prob - ',current_image);
            p = p + 1;
            
            figure(p)
            h = heatmap(Prob_Map(:,:,3),'ColorLimits',[0 1]);
            h.Title = strcat('k=5 - Cluster 3 Prob - ',current_image);
            p = p + 1;
            
            figure(p)
            h = heatmap(Prob_Map(:,:,4),'ColorLimits',[0 1]);
            h.Title = strcat('k=5 - Cluster 4 Prob - ',current_image);
            p = p + 1;
            
            figure(p)
            h = heatmap(Prob_Map(:,:,5),'ColorLimits',[0 1]);
            h.Title = strcat('k=5 - Cluster 5 Prob',current_image);
            p = p + 1;
        
        end
    end
    
    % Analysis of inputs for debugging
    figure(p)
    subplot(2, 2, 1);
    imshow(Image_A);
    title(strcat('Original Color Image - ',current_image));
    subplot(2, 2, 2);
    imshow(red);
    title('Red Image');
    subplot(2, 2, 3);
    imshow(green);
    title('Green Image');
    subplot(2, 2, 4);
    imshow(blue);
    title('Blue Image');
    p = p + 1;
    
    % Plotting of clusters and results
    % figure(p)
    % subplot(2, 2, 1);
    % imshow(Image_A);
    % title(strcat('Original Image - ',current_image));
    % subplot(2, 2, 2);
    % image(Plot_Map_2);
    % title('2 Clusters');
    % subplot(2, 2, 3);
    % image(Plot_Map_3);
    % title('3 Clusters');
    % subplot(2, 2, 4);
    % image(Plot_Map_5);
    % title('5 Clusters');
    % p = p + 1;

end