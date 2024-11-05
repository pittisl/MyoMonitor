% This is the main algorithm
% This has to use the data collected by MyoTester App

maxChanTap = 4;
ChanTapSamples = 300;

filepath = uigetdir('F:\Matlab_Files_2019b\muscle_v3\GroupData','Select Folder for PCM files');
if filepath == 0
    return
end
listing = dir([filepath,'\3Hold\*.pcm']);

% featureMatrix = zeros(length(listing), (maxChanTap - 1) * ChanTapSamples * 1); %
featureMatrix = zeros(length(listing), 11); %
% featureMatrix = zeros(length(listing), ChanTapSamples); %
% acousticData = zeros(3,ChanTapSamples,length(listing));
for fileIndex = 1:length(listing)
    %% Read each data file
    pcm_data = ReadAudioFile([filepath, '\3Hold\', listing(fileIndex).name]);
    rawFeatures = preprocessing(pcm_data, maxChanTap);
    
    %% Feature extraction
    readyFeatures = featureExtraction(rawFeatures, ChanTapSamples);
%     acousticData(:,:,fileIndex) = readyFeatures;
    featureMatrix(fileIndex, :) = readyFeatures; %% Rows are observations
end

% featureMatrix = mapminmax(featureMatrix', 0, 1);
% featureMatrix = featureMatrix';
% figure,plot(featureMatrix)

%% Classification
% p2=fatigueClassify(featureMatrix);

labels = dbscan(kD, (kD(1,2)+kD(1,3)+kD(2,3))./12, 3, 'Distance', 'precomputed');
featureMatrix(find(labels == -1), :) = [];
kD = pdist2(featureMatrix,featureMatrix,'euc');
labels = dbscan(kD, (kD(1,2)+kD(1,3)+kD(2,3))./30, 3, 'Distance', 'precomputed');