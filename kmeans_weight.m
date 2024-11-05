% This is to calculate k-means clustering for single muscle workout
% Initialization  DO THIS!
k_means_tmp = [];
% Get k means observations one by one through the data
data = bpsk_0227_10_09d0p;       % CHANGE THIS TO VARIABLE NAME
k_means_tmp_ = k_means_tmp;
k_means_tmp = kmeans_segformat(data, k_means_tmp); % fill in the variable name
k_means_tmp = k_means_tmp_;

% calculate kmeans
k_means_tmp = k_means_tmp';
% option: initialize with two mean centroids
cluster1 = mean(k_means_tmp(1:5,:));
cluster2 = mean(k_means_tmp(21:25,:));
[k_idx, k_Cent] = kmeans(k_means_tmp, 2, 'Start', [cluster1;cluster2]);
% ground truth clustering vector
k_idx_ground = [repelem(1,20),repelem(2,20)];
% Accuracy
cnt = 0;
for i = 1:30
    if k_idx(i) == k_idx_ground(i)
        cnt = cnt + 1;
    end
end
cnt./30     % print accuracy

% calculate knn
% take out training set and test set
knn_train = [k_means_tmp(1:20,:);k_means_tmp(21:45,:);k_means_tmp(51:55,:)];
knn_test = [k_means_tmp(46:50,:);k_means_tmp(56:70,:);k_means_tmp(71:80,:)];
knn_train_label = [repelem(1,20),repelem(2,25),repelem(1,5)];
knn_test_label = [repelem(2,5),repelem(1,15),repelem(2,10)];
knnMd = fitcknn(knn_train, knn_train_label);
% test 
label = predict(knnMd,knn_test);
% Accuracy
cnt = 0;
for i = 1:length(label)
    if label(i) == knn_test_label(i)
        cnt = cnt + 1;
    end
end
cnt./30     % print accuracy