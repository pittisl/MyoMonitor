% This is the main algorithm
% This has to use the data collected by MyoTester App
% Classify non-fatigue vs. fatigue
% Data are all the beginning 5 vs. last 5

ChanTapSamples = 300;

%% Non-fatigue
filepath = uigetdir('F:\Matlab_Files\muscle_v3\GroupData','Select Folder for non-fatigue PCM files');
if filepath == 0
    return
end
listing = dir([filepath,'\*.pcm']);

featureMatrix_nonFatigue = zeros(length(listing), ChanTapSamples);
label_nonFatigue = zeros(length(listing), 1);
for fileIndex = 1:length(listing)
    %% Read each data file
    pcm_data = ReadAudioFile([filepath, '\', listing(fileIndex).name]);
    rawFeatures = preprocessing(pcm_data, 2);
    
    %% Feature extraction
    readyFeatures = featureExtraction(rawFeatures, ChanTapSamples);
    featureMatrix_nonFatigue(fileIndex, :) = readyFeatures; %% Rows are observations
end

%% Fatigue
filepath = uigetdir('F:\Matlab_Files\muscle_v3\GroupData','Select Folder for fatigue PCM files');
if filepath == 0
    return
end
listing = dir([filepath,'\*.pcm']);

featureMatrix_Fatigue = zeros(length(listing), ChanTapSamples);
label_Fatigue = ones(length(listing), 1);
for fileIndex = 1:length(listing)
    %% Read each data file
    pcm_data = ReadAudioFile([filepath, '\', listing(fileIndex).name]);
    rawFeatures = preprocessing(pcm_data, 2);
    
    %% Feature extraction
    readyFeatures = featureExtraction(rawFeatures, ChanTapSamples);
    featureMatrix_Fatigue(fileIndex, :) = readyFeatures; %% Rows are observations
end

% Entire feature matrix
rawfeatures = [featureMatrix_Fatigue; featureMatrix_nonFatigue];
labels = [label_Fatigue; label_nonFatigue];

%% Feature extraction
% normalization
for i = 1:size(rawfeatures, 1)
    rawfeatures(i, :) = (rawfeatures(i, :) - mean(rawfeatures(i, :)))/var(rawfeatures(i, :));
end

% time domain features
featuresTime = zeros(size(rawfeatures, 1), 6);
for i = 1:size(rawfeatures, 1)
    featuresTime(i, :) = [mean(rawfeatures(i,:)), ...
                            var(rawfeatures(i,:)), ...
                            skewness(rawfeatures(i,:)), ...
                            min(rawfeatures(i,:)), ...
                            max(rawfeatures(i,:)), ...
                            mean(abs(diff(sign(rawfeatures(i,:)))))];
end

% frequency domain features
featuresSpec = zeros(size(rawfeatures));
for i = 1:size(rawfeatures, 1)
    featuresSpec(i, :) = fft(rawfeatures(i, :));
end

featuresFreq = zeros(size(featuresSpec, 1), 5);
for i = 1:size(featuresSpec, 1)
    featuresFreq(i, :) = [mean(abs(featuresSpec(i,:))), ...
                            var(abs(featuresSpec(i,:))), ...
                            skewness(abs(featuresSpec(i,:))), ...
                            min(abs(featuresSpec(i,:))), ...
                            max(abs(featuresSpec(i,:))), ...
                            sum(featuresSpec(i,:).*conj(featuresSpec(i,:)))];
end

features = featuresFreq;

%% Classify
correct = 0;
for validationCnt = 1:100
    % pre-permutation
    randomRows = randperm(70);
    randomRows = randomRows(1:50);
    trainSet = features(randomRows,:);
    testSet = features;
    testSet(randomRows,:) = [];
    trainLabels = labels(randomRows);
    testLabels = labels;
    testLabels(randomRows) = [];
    % 1.kNN
%     Mdl = fitcknn(trainSet, trainLabels);
%     label_pre = predict(Mdl, testSet);

    % 2.SVM
%     Md1 = fitcsvm(trainSet, trainLabels, 'KernelFunction', 'rbf');
%     label_pre = predict(Md1, testSet);
    
    % 3.NB
%     Md1 = fitcnb(trainSet, trainLabels);
%     label_pre = predict(Md1, testSet);
    
    % 4.Decision tree
%     Md1 = fitctree(trainSet, trainLabels);
%     label_pre = predict(Md1, testSet);

    % 5.Discriminant Analysis
    Md1 = fitcdiscr(trainSet, trainLabels, 'DiscrimType', 'pseudoLinear');
    label_pre = predict(Md1, testSet);

    % Accuracy
    for i = 1:length(label_pre)
        if testLabels(i) == label_pre(i)
            correct = correct + 1;
        end
    end
end
accuracy = correct / (length(label_pre)*100);
sprintf('Accuracy is %2.2f%%', accuracy*100)