% This is to classify fatigue
% This has to use the data collected by MyoTester App
function pFeature = fatigueClassify(features)
    pFeature = 1;
    
    %% Difference between every two data
    
%     featureDiff = zeros(size(features, 1) - 1, size(features, 2));
%     for i = 1:size(featureDiff)
%         featureDiff(i, :) = sqrt(sum((features(i + 1, :) - features(i, :)) .^ 2));
%     end
% %     figure, 
%     plot(sum(featureDiff, 2))
    %% remove some outliers
    %% 1. MAD
%     featureMean = mean(features);
%     featureMeanDist = zeros(size(features, 1), 1);
%     for i = 1:size(features, 1)
%         featureMeanDist(i) = sqrt(sum((features(i, :) - featureMean) .^ 2));
%     end
%     outlierIndex = isoutlier(featureMeanDist, 'movmedian', 5);
%     features(outlierIndex, :) = [];
%     figure, plot(featureMeanDist)
%     disp(['Outliers: ', num2str(sum(outlierIndex))])
    
    %% smooth
%     for i = 1:size(features, 2)
%        features(:, i) = smooth(features(:, i), 5);
%     end
    
    %% Linearization
%     num = 1:size(features, 1);
%     pFeature = zeros(1, size(features, 2));
%     for i = 1:size(features, 2)
%         [p, S] = polyfit(num', features(:, i), 1);
%         pFeature(i) = p(1);
%         [features(:, i), ~] = polyval(p, num, S);
%     end
%     disp(['Mean slope: ', num2str(mean(pFeature))])
    
    %% kmeans
%     numOfObservations = size(features, 1);
%     feedIndices = zeros(1, numOfObservations);
%     for i = 1:floor(numOfObservations/2)
%         feedIndices(2*i - 1) = i;
%         feedIndices(2*i) = numOfObservations - i + 1;
%     end
%     if mod(numOfObservations, 2) == 1
%         feedIndices(end) = ceil(numOfObservations/2);
%     end
%     
%     featuresKmeans = features(feedIndices, :);
    %% trials
    % 1.
%     startKmeans = [features(1,:),features(end,:), ...
%         features(2,:), features(end-1,:), ...
%         features(3,:), features(end-2,:)];
    % 2.
%     startKmeans(:,:,1) = [features(1,:);features(end,:)];
%     startKmeans(:,:,2) = [features(2,:);features(end-1,:)];
%     startKmeans(:,:,3) = [features(3,:);features(end-2,:)];
    
    % 3.
%     startKmeans = [features(1,:); features(end,:); ...
%                     features(2,:); features(end-1,:); ...
%                     features(3,:); features(end-2,:)];
%     
%     [idx_, startC] = kmeans(startKmeans,2, 'Replicates', 50);
%     figure, plot(idx_)

    % 4. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     startKmeans = [mean(features(1:5, :)); ...
% %         mean(features(round(size(features, 1)/2)-1 : round(size(features, 1)/2) + 1, :)); ...
%         mean(features(end - 4:end, :)) ];
%     for i = 1:50
%         startC(:,:,i) = startKmeans; 
%     end
%     [idx, C] = kmeans(features(4: end-3, :), 2, 'Replicates', 50, 'Start', startC);
% %     [idx, C] = kmeans(features, 2, 'Replicates', 50);
%     if idx(1) == 2 % make the beginning always the non-fatigue
%         for i = 1:length(idx)
%             if (idx(i) == 2)
%                 idx(i) = 1;
%             else
%                 idx(i) = 2;
%             end
%         end
%     end
%     idx = [1;1;1;idx;2;2;2];
%     
%     %% plot
%     figure, plot(idx, 'xr', 'MarkerSize', 15)
%     if idx(1) == 1
%         names = {'Non-fatigue'; ...
% %             'Transition'; ...
%                 'Fatigue'};
%     else
%         names = {'fatigue'; ...
% %             'Transition'; ...
%                 'Non-fatigue'};
%     end
%     set(gca,'ytick',1:2,'yticklabel',names)
%     xlabel('Workout index')
%     title('K-means cluster results')
%     grid on
%     axis([-inf inf 0 3]);
    
%% DBSCAN


end