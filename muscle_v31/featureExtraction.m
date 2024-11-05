% This is to extract features from raw channel estimation
% This has to use the data collected by MyoTester App
function readyFeatures = featureExtraction(rawFeatures, chanTapSamples)
    %% Our sampling rate is 100Hz, so interpolate the data to 300 points
    interpFeatures = zeros(size(rawFeatures, 1), chanTapSamples);
    for i = 1:size(rawFeatures, 1)
        interpFeatures(i, :) = interp1(1:size(rawFeatures, 2), rawFeatures(i, :), ...
            linspace(1, size(rawFeatures, 2), chanTapSamples));
    end
    

    %% Normalization
    for i = 1:size(interpFeatures, 1)
%         interpFeatures(i, :) = (interpFeatures(i, :) - mean(interpFeatures(i, :)))/var(interpFeatures(i, :));
%         interpFeatures(i, :) = mapminmax(interpFeatures(i, :));
    end
    
    
    %% Spectrum
    specFeatures = zeros(size(interpFeatures));
    for i = 1:size(specFeatures, 1)
        specFeatures(i, :) = abs(fft(interpFeatures(i, :)));
    end
    
    %% time domain features
    timeRawFeatures = interpFeatures(1,:);
%     featuresTime = zeros(1, 6);
%     for i = 1:6
        featuresTime = [mean(timeRawFeatures), ...
                            var(timeRawFeatures), ...
                            skewness(timeRawFeatures), ...
                            min(timeRawFeatures), ...
                            max(timeRawFeatures)];
%     end

    %% frequency domain features
    featuresSpec = fft(timeRawFeatures);

%     featuresFreq = zeros(size(featuresSpec, 1), 5);
%     for i = 1:size(featuresSpec, 1)
        featuresFreq = [mean(abs(featuresSpec)), ...
                                var(abs(featuresSpec)), ...
                                skewness(abs(featuresSpec)), ...
                                min(abs(featuresSpec)), ...
                                max(abs(featuresSpec)), ...
                                sum(featuresSpec.*conj(featuresSpec))];
%     end

    
    %% Return ready features
%     readyFeatures = reshape(interpFeatures', 1, []);
%     readyFeatures = interpFeatures(1,:);
%     readyFeatures = mean(interpFeatures);
%     figure, plot(readyFeatures)
    readyFeatures = [featuresTime, featuresFreq];
    
    %% Plot
%     figure
%     for i = 1:size(specFeatures, 1)
%         subplot(size(specFeatures, 1),1,i), hold on
%         plot(specFeatures(i,:), 'LineWidth', 1)
%     end
%     figure, plot(readyFeatures)
end