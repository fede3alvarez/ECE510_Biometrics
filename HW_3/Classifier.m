close all
clear all;
clc;

% Code will be available at
% https://github.com/fede3alvarez/ECE510_Biometrics

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Parameters of interest                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k_values = 3; %[2 3 5];   % Array with the # of k-clusters to test
k_repeats = 100;          % # of times to repeat kmeans
prob_lim = [0, 0.2, 0.6]; %curr_prob_lim

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
                    %'Hand01.jpg'
                    %'Hand02.jpg'
                    %'Hand03.jpg'
                    'Hand04.jpg'
                    ];

for m = 1:size(Available_images,1)
    current_image = Available_images(m,:);
    
    % Note: jpeg images are read directly as RGB  matrices
    %    https://www.mathworks.com/matlabcentral/answers/23119-jpg-to-rgb-image
    [Image_A, Map_A] = imread(current_image);

    % Convert to Grayscale
    Gray_Map_A = rgb2gray(Image_A);

    % Massage data
    data = double([Gray_Map_A]);
    [rows, cols] = size(data);
    data_list = data(:);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %              K-Means Clustering                %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % In this section we setup the numbers of K cluster to be used
    
    for k = k_values
        % From: Homework Prompt:
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
                init_cond = [init_cond; randsample(data_list,1)];
            end        
            
            % Sort clusters based on Darker (k=1) TO Lighter (k=3)
            % Note: In Binary Image and Grayscale Black=0, White=1
            init_cond = sortrows(init_cond,'descend');
        
            % Run k-means
            [idx, C] = kmeans(data_list, k,'start', init_cond);
        
            
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
        
        end  % End of k_repeats
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                 Find Boundary                  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

        BW = imbinarize(Gray_Map_A);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                  Plots &                       %
        %                SPECIFICATIONS                  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Probabilities Threesholds
        prob_lim = [0, 0.2, 0.6]; %curr_prob_lim

        for lim = prob_lim

            % Remove any probability below ther probability limnit
            Prob_Map(Prob_Map<lim) = 0;

            % Plot image and boundaries
            figure(p)
            subplot(2, 3, 1);
            imshow(Image_A);
            title(strcat('Original Color Image - ',current_image));
    
            subplot(2, 3, 2);
            imshow(BW);
            title(strcat('Binary Image - ',current_image));
    

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %              SPECIFICATIONS               %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % Assumptions
            % Since the k-means has been organized from darkest to lighter:
            %   k=1 is assumed to be black / background
            %   k=3 is the lighter / hand
            %   k=2 is the mix / boundaries
            %
            % Since the binary image is used as the ground-truth
            %   Black = 0 = FALSE
            %   White = 1 = TRUE

            max_prob = max(Prob_Map,[],3);

            % Plot Values
            % True positive is where:
            %   Binary White (BW=1) matches Result White (K=3)
            idx = Prob_Map(:,:,3);
            idx(idx == max_prob) = 1; 
            temp = (idx == BW);
            TP = sum(temp(:));

            % True negative is where:
            %   Binary Black (BW=0) matches Result Black (K=1)
            idx = Prob_Map(:,:,1);
            idx(idx == max_prob) = 1; 
            temp = (idx == ~BW);
            TN = sum(temp(:));

            % False positive is where:
            %    Binary White (BW=1) matches Result Black (K=1)
            idx = Prob_Map(:,:,1);
            idx(idx == max_prob) = 1;
            temp = (idx == BW);
            FP = sum(temp(:));

            % False negative is where:
            %   Binary Black (BW=0) matches Result White (K=3)
            idx = Prob_Map(:,:,3);
            idx(idx == max_prob) = 1; 
            temp = (idx == ~BW);
            FN = sum(temp(:));
            
            % Total count of Binary White
            temp = (BW(:) == 1);
            POS = sum(BW(:) == 1);
            
            % Total count of Binary Black
            NEG = sum(BW(:) == 0);

            ACC = (TP+TN) / (POS+NEG);  % ACCURACY
            TPR = TP / (TP+FN);         % SENSITIVITY (TRUE POS RATE)
            FPR = FP / (FP+TN);         % FALSE POS RATE
            TNR = TN / (FP+TN);         % SPECIFICITY

            % Plotting information
            subplot(2, 3, 3);
            n = newline;

            Explanation = strcat('Accuracy: ', sprintf('%0.3f',ACC),'\n',...
                                 'TP Rate: ', sprintf('%0.3f',TPR),'\n',...
                                 'FP Rate: ', sprintf('%0.3f',FPR),'\n',...
                                 'TN Rate:', sprintf('%0.3f',TNR));
            Explanation = compose(Explanation);
            text(0,0.5,Explanation); axis off


            subplot(2, 3, 4);
            h1 = heatmap(Prob_Map(:,:,1));
            h1.title(strcat('Prob Lim:',...
                            sprintf('%0.3f',lim),...
                            ' Cluster:1-Black/Background'));

            subplot(2, 3, 5);
            h1 = heatmap(Prob_Map(:,:,2));
            h1.title(strcat('Prob Lim:',...
                            sprintf('%0.3f',lim),... 
                            ' Cluster:2-Boundary'));
    
            subplot(2, 3, 6);
            h1 = heatmap(Prob_Map(:,:,1));
            h1.title(strcat('Prob Lim:',...
                             sprintf('%0.3f',lim),...
                             ' Cluster:3-White/Hand'));
    


            p = p + 1; 
            % % Plotting of Probability Map
            % figure(p)
            % 
            % for j = 1:k 
            % h = heatmap(Prob_Map(:,:,1));
            % h.Title = strcat('Cluster:',int2str(j),...
            %                  ' Prob Lim:',int2str(lim),...
            %                  ' Image:',current_image);
            % p = p + 1;
            %end
        end
    end % End of Images Loop
    


end