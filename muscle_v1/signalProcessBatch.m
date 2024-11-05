% This is to process the received data IN FOLDER
% clear;clc;
function fatigueValues = signalProcessBatch()
% Parameters
seqSymbols = [-1,-1, 1,-1,-1, ...
               1,-1, 1, 1, 1, ...
              -1,-1,-1,-1, 1, ...
              -1,-1,-1, 1,-1, ...
              -1, 1,-1, 1, 1, 1];
fs = 48000;
fc = 20000;
c = 340;
codeUpsampleRate = 4;
P = 16;
L = 10;
frameLength = 50;
% thresholdStartPoint = 0.6917;
guardSamples = 10000;
perTapDist = 100 * c/(2 * fc/codeUpsampleRate); % in centimeter
load('M_mat.mat'); % load matrix M 
colorMap = [0 0 0;
            0 0 1;
            0 1 0;
            0.5 0 0.5;
            1 0 0;
            1 0 1;
            0.5 0.5 0;
            0 0 0.5;
            0 0.5 0;
            0 0.5 0.5;
            0.5 0 0];
dura_stat = 3; dura_power = 0; dura_lift = 3; dura_hold = 1; dura_put = 3; % standard test timing
dura_total = dura_stat + dura_power + dura_lift + dura_hold + dura_put;
        
%% Read data from PCM file IN FOLDER
filepath = uigetdir('F:\Matlab_Files\muscle_v1\Group_Data','Select Folder for PCM files');
if filepath == 0
    return
end
listing = dir([filepath,'\*.pcm']);
repreValueFatigue = zeros(1,length(listing));  %%%% COMMENT FOR 1 feature
% ----
% repreValueFatigue = zeros(2,length(listing));
% ----
sumChanEstiValueReal = zeros(L+1,length(listing));
sumChanEstiValueImag = zeros(L+1,length(listing));
for fileIndex = 1:length(listing)
%     pcm_data = ReadAudioFile([filePath, listing(fileIndex).name]); % ReadAudioFile should be under same directory
    pcm_data = ReadAudioFile([filepath, '\', listing(fileIndex).name]);
    % figure, plot(pcm_data)
    load('ref_code.mat'); % load reference code file in variable 

    %% Compute autocorrelation to detect the start point
    corr = zeros(1, guardSamples);
    for i = 1:guardSamples
        R = corrcoef(pcm_data(i:i+length(y)-1),y);
        corr(i) = R(2,1);
    end
    % figure,plot(abs(corr))

    %% Cut the data from the start point
    [~, startPoint] = max(abs(corr));
    pcm_data = pcm_data(startPoint:end)';

    %% Down-conversion - Multiplication
    t = 0:1/fs:length(pcm_data)/fs;
    real_cos = sqrt(2) * cos(2*pi*fc*t);
    imag_sin = -sqrt(2) * sin(2*pi*fc*t);
    real_base = real_cos(1:end-1) .* pcm_data;
    imag_base = imag_sin(1:end-1) .* pcm_data;
    % real_base = real_base./abs(real_base);   % normalize
    % imag_base = imag_base./abs(imag_base);

    %% Taking average over cycle of sample rate
    samplesPerCode = fs*codeUpsampleRate/fc;
    real_mean = zeros(1, floor(length(real_base)/samplesPerCode));
    imag_mean = zeros(1, floor(length(imag_base)/samplesPerCode));
    for i = 1 : length(real_base)/samplesPerCode
        real_mean(i) = mean(real_base(round((i-1)*samplesPerCode + 1):round(i * samplesPerCode)));
        imag_mean(i) = mean(imag_base(round((i-1)*samplesPerCode + 1):round(i * samplesPerCode)));
    end

    %% Complex baseband signal
    r_complex = real_mean + 1j .* imag_mean;

    %% Take out yL+1 to yL+P and form a vector yL
    r_complex = r_complex(1:end - mod(length(r_complex), frameLength));
    y_L_matrix = reshape(r_complex, 50, []);
    y_L_matrix = y_L_matrix(L+1: L+P, :);       % each column is an observation

    %% Compute channel estimation h^
    h_ = (M' * M) \ M' * y_L_matrix;
    [l_r1,l_c1] = size(h_);

    % plot the complex value of h_
%     figure, hold on
%     for i = 1 : L+1
% %         subplot(4,3,i), 
%         hold on
%             plot(h_(i, :), ...
%             'Color', colorMap(i,:), 'LineWidth', 2)
%         title([num2str(perTapDist * (i-1)), ' to ',num2str(perTapDist * i), 'cm'])
% %         legend('show', 'Location', 'southeast')
%     end
%     suptitle('IQ traces from different channel tap')
%     legend('show', 'Location', 'southeast')

    %% Compute the difference h^t - h^t-1
    h_d_ = zeros(l_r1,l_c1-1);
    for i = 1 : l_c1-1
        h_d_(:,i) = h_(:,i + 1) - h_(:,i);
    end
%     l_c1 = l_c1 - 1;?% if use h_d_ later
    % plot the complex value of h_d_
    % figure, hold on
    % for i = 1 : L+1
    %     subplot(4,3,i),plot(h_d_(i, :), 'Color', colorMap(i,:), ...
    %         'DisplayName', [num2str(perTapDist * (i-1)), ' to ',num2str(perTapDist * i), 'cm'], ...
    %         'LineWidth', 2)
    %     legend('show', 'Location', 'northwest')
    % end
    % suptitle('IQ traces from different channel tap (after static removal)')
    % legend('show', 'Location', 'northwest')

    %% Compute the summation of real and imaginary components
%     sumChanEstiValueReal(:,fileIndex) = sum(abs(h_'));
%     sumChanEstiValueImag(:,fileIndex) = sum(imag(h_'));

    %% Filter out frequency components outside 8-12Hz
%     h_spec = fft(h_');  % the transpose'
%     % Take out spectrum of 8-12 Hz
%     sampleRate = fc/codeUpsampleRate/frameLength;
%     spec_start = l_c1 * (8/sampleRate);
%     spec_end = l_c1 * (12/sampleRate);
%     % filter out other frequency
%     h_spec(1:floor(spec_start - 1), :) = 0;
%     h_spec(ceil(spec_end + 1): end, :) = 0;
%     h_ = ifft(h_spec);
%     h_ = h_';

    %% plot each tap in different workout segments
    period_stat = 1 : round(l_c1*dura_stat/dura_total);     % segments
    period_power = round(l_c1*dura_stat/dura_total) + 1 : ...
        round(l_c1*(dura_stat+dura_power)/dura_total);
    period_lift = round(l_c1*(dura_stat+dura_power)/dura_total) + 1 : ...
        round(l_c1*(dura_stat+dura_power+dura_lift)/dura_total);
    period_hold = round(l_c1*(dura_stat+dura_power+dura_lift)/dura_total) + 1 : ...
        round(l_c1*(dura_stat+dura_power+dura_lift + dura_hold)/dura_total);
    period_put = round(l_c1*(dura_stat+dura_power+dura_lift + dura_hold)/dura_total) + 1 : ...
        round(l_c1*(dura_stat+dura_power+dura_lift + dura_hold + dura_put)/dura_total);
    % plot
    centerIcom = zeros(L+1, 5);
    centerQcom = zeros(L+1, 5);
    sumDistance = zeros(L+1, 5); %%%% COMMENT FOR 1 feature
    % ----
%     sumDistance = zeros(L+1, 10);  
    % ----
%     figure,
    for i = 1 : L+1
%         subplot(4,3,i), hold on
%         plot(h_(i, period_stat), ...
%             'LineWidth', 1, 'DisplayName', 'static')
%         plot(h_(i, period_power), ...
%             'LineWidth', 1, 'DisplayName', 'generate')
%         plot(h_(i, period_lift'), ...
%             'LineWidth', 1, 'DisplayName', 'lift-up')
%         plot(h_(i, period_hold), ...
%             'LineWidth', 1, 'DisplayName', 'hold')
%     %     plot(h_(i, period_put), ...
%     %         'LineWidth', 1, 'DisplayName', 'lift-down')
%         title([num2str(perTapDist * (i-1)), ' to ',num2str(perTapDist * i), 'cm'])
%         legend('show', 'Location', 'northeast')
% %         calculate the centroid of the point set in each workout section and
% %         sum the distance of all the points
        centerIcom(i,:) = [mean(real(h_(i, period_stat))), mean(real(h_(i, period_power))), ...
            mean(real(h_(i, period_lift))), mean(real(h_(i, period_hold))), mean(real(h_(i, period_put)))];
        centerQcom(i,:) = [mean(imag(h_(i, period_stat))), mean(imag(h_(i, period_power))), ...
            mean(imag(h_(i, period_lift))), mean(imag(h_(i, period_hold))), mean(real(h_(i, period_put)))];
        sumDistance(i,:) = [sumDistanceFunc(h_(i, period_stat), centerIcom(i,1), centerQcom(i,1)), ...
            sumDistanceFunc(h_(i, period_power), centerIcom(i,2), centerQcom(i,2)), ...
            sumDistanceFunc(h_(i, period_lift), centerIcom(i,3), centerQcom(i,3)), ...
            sumDistanceFunc(h_(i, period_hold), centerIcom(i,4), centerQcom(i,4)), ...
            sumDistanceFunc(h_(i, period_put), centerIcom(i,5), centerQcom(i,5))];    %%%% COMMENT FOR 1 feature
        % ----
%         [sumDistanceTemp_stat, sumDistanceTemp2_stat]= sumDistanceFunc(h_(i, period_stat), centerIcom(i,1), centerQcom(i,1));
%         [sumDistanceTemp_power, sumDistanceTemp2_power]= sumDistanceFunc(h_(i, period_power), centerIcom(i,2), centerQcom(i,2));
%         [sumDistanceTemp_lift, sumDistanceTemp2_lift]= sumDistanceFunc(h_(i, period_lift), centerIcom(i,3), centerQcom(i,3));
%         [sumDistanceTemp_hold, sumDistanceTemp2_hold]= sumDistanceFunc(h_(i, period_hold), centerIcom(i,4), centerQcom(i,4));
%         [sumDistanceTemp_put, sumDistanceTemp2_put]= sumDistanceFunc(h_(i, period_put), centerIcom(i,5), centerQcom(i,5));
%         sumDistance(i,:) = [sumDistanceTemp_stat, sumDistanceTemp_power, sumDistanceTemp_lift, ...
%                             sumDistanceTemp_hold, sumDistanceTemp_put, sumDistanceTemp2_stat, ...
%                             sumDistanceTemp2_power, sumDistanceTemp2_lift, sumDistanceTemp2_hold, ...
%                             sumDistanceTemp2_put];
        % ----
    end
    % suptitle('IQ traces from different channel tap with workout label')
    
    %% compute the representative values (average 1-3 minus average 4-11)
%     repreValueFatigue(1, fileIndex) = mean(sumDistance(1:3, 4)) - mean(sumDistance(4:L+1, 4));
    % ----
%     repreValueFatigue(2, fileIndex) = mean(sumDistance(1:3, 9));% - mean(sumDistance(4:L+1, 9));
    % ----
    repreValueFatigue(fileIndex) =mean(sumDistance(1:3, 4)) - mean(sumDistance(4:L+1, 4)); %% COMMENT FOR 1 feature
    
    %% plot the degree of concentration
%     figure,hold on
%     for i = 1 : L+1
%         plot(sumDistance(i,:), 'Color', colorMap(i,:), ...
%             'DisplayName', [num2str(perTapDist * (i-1)), ' to ',num2str(perTapDist * i), 'cm'], ...
%             'LineWidth', 2)
%         legend('show', 'Location', 'northeast')
%     end
%     title('IQ traces level of concentration at each workout stage')

    %% Plot distance derived from phase of channel estimation
    % h_abs = zeros(l_r1, l_c1);
    % figure, hold on
    % for i = 1 : L+1
    %     plot(linspace(0, l_c1 * frameLength * samplesPerCode / fs, l_c1)+startPoint/fs, ...
    %         (smooth(phase(h_(i, :)))/(2*pi))*(c/fc)*0.5*100, ...
    %         'Color', colorMap(i,:), 'DisplayName', [num2str(perTapDist * (i-1)), ' to ',num2str(perTapDist * i), 'cm'], ...
    %         'LineWidth', 2)
    %     h_abs(i, :) = (smooth(phase(h_(i, :)))/(2*pi))*(c/fc)*0.5*100;
    % end
    % title('Distance Channel Tap')
    % xlabel('Time(s)'), ylabel('Distance(cm)')
    % legend('show', 'Location', 'southwest')

    %% Plot raw phase of channel estimation
%     figure, hold on
%     for i = 1 : L+1
%         subplot(4,3,i), plot(linspace(0, l_c1 * frameLength * samplesPerCode / fs, l_c1)+startPoint/fs, ...
%             smooth(angle(h_(i, :))), ...
%             'Color', colorMap(i,:), 'DisplayName', [num2str(perTapDist * (i-1)), ' to ',num2str(perTapDist * i), 'cm'], ...
%             'LineWidth', 2)
%         xlabel('Time(s)')
%         legend('show', 'Location', 'southwest')
%     end
%     suptitle('Raw Phase Channel Tap')
%     xlabel('Time(s)')
%     legend('show', 'Location', 'southwest')

    %% Plot amplitude of channel estimation
    % figure, hold on
    % for i = 1 : L+1
    %     subplot(4,3,i), plot(linspace(0, l_c1 * frameLength * samplesPerCode / fs, l_c1)+startPoint/fs, ...
    %         smooth(abs(h_(i, :))), ...
    %         'Color', colorMap(i,:), 'DisplayName', [num2str(perTapDist * (i-1)), ' to ',num2str(perTapDist * i), 'cm'], ...
    %         'LineWidth', 2)
    %     xlabel('Time(s)')
    %     legend('show', 'Location', 'southwest')
    % end
    % suptitle('Raw Amplitude Channel Tap')
    % xlabel('Time(s)')
    % legend('show', 'Location', 'northeast')

    %% Try MTI method from Yuqi
    % for i = 1: L+1
    %     [h_(i,:), ~] = ModifiedMTI(h_(i,:), 4);
    % end
    % h_d_ = h_; l_c1 = l_c1+1;

    % figure, plot(h_d_(1,100:200))
    % plot every channel tap
    % for i = 1 : L+1
    %     figure, plot(linspace(0, (l_c1-1) * frameLength * samplesPerCode / fs, l_c1-1)+startPoint/fs, (smooth(phase(h_d_(i, :)))/(2*pi))*(c/fc)*0.5)
    %     title(['Channel Tap ', num2str(i)])
    %     axis([-inf inf -0.2 0.2]);
    % end

    %% Plot distance derived from phase of channel estimation difference
    % figure, hold on
    % for i = 1 : L+1
    %     plot(linspace(0, (l_c1-1) * frameLength * samplesPerCode / fs, l_c1-1)+startPoint/fs, ...
    %         (smooth(phase(h_d_(i, :)))/(2*pi))*(c/fc)*0.5*100, ...
    %         'Color', colorMap(i,:), 'DisplayName', [num2str(perTapDist * (i-1)), ' to ',num2str(perTapDist * i), 'cm'], ...
    %         'LineWidth', 2)
    %     xlabel('Time(s)'), ylabel('Distance(cm)')
    %     legend('show', 'Location', 'southwest')
    % end
    % suptitle('Distance from Channel Tap (after static objects removal)')
    % xlabel('Time(s)'), ylabel('Distance(cm)')
    % legend('show', 'Location', 'southwest')

    %% Plot Phase of channel estimation difference
    % figure, hold on
    % for i = 1 : L+1
    %     subplot(4,3,i), plot(linspace(0, (l_c1-1) * frameLength * samplesPerCode / fs, l_c1-1)+startPoint/fs, ...
    %         smooth(angle(h_d_(i, :))), ...
    %         'Color', colorMap(i,:), 'DisplayName', [num2str(perTapDist * (i-1)), ' to ',num2str(perTapDist * i), 'cm'], ...
    %         'LineWidth', 2)
    %     xlabel('Time(s)')
    %     legend('show', 'Location', 'southwest')
    % end
    % suptitle('Phase from Channel Tap (after static objects removal)')
    % xlabel('Time(s)')
    % legend('show', 'Location', 'southwest')

    % %% Plot amplitude of channel estimation difference
    % figure,
    % for i = 1 : L+1
    %     subplot(4,3,i), plot(linspace(0, (l_c1-1) * frameLength * samplesPerCode / fs, l_c1-1)+startPoint/fs, ...
    %         smooth(abs(h_d_(i, :)), 20), 'Color', colorMap(i,:), ...
    %         'DisplayName', [num2str(perTapDist * (i-1)), ' to ',num2str(perTapDist * i), 'cm'])
    %     xlabel('Time(s)')
    %     legend('show', 'Location', 'northeast')
    % end
    % suptitle('Amplitude from Channel Tap (after static objects removal)')

    %% PCA
    % h_abs = (smooth(phase(h_))/(2*pi))*(c/fc)*0.5*100;
    % h_n = mapminmax(h_abs);
    % figure ,plot(h_n')
    % % figure,
    % % for i = 1 : L+1
    % %     subplot(4,3,i), plot(linspace(0, (l_c1) * frameLength * samplesPerCode / fs, l_c1)+startPoint/fs, ...
    % %         h_n(i,:), 'Color', colorMap(i,:), ...
    % %         'DisplayName', [num2str(perTapDist * (i-1)), ' to ',num2str(perTapDist * i), 'cm'])
    % %     xlabel('Time(s)')
    % %     legend('show', 'Location', 'northeast')
    % % end
    % [coeff,score,latent] = pca(h_n');
    % % pca_com = score(:,4);
    % figure, hold on
    % for i = 1 : L+1
    %     plot(linspace(0, (l_c1) * frameLength * samplesPerCode / fs, l_c1)+startPoint/fs, ...
    %         score(:,i), 'Color', colorMap(i,:), ...
    %         'DisplayName', [num2str(perTapDist * (i-1)), ' to ',num2str(perTapDist * i), 'cm'])
    % %     xlabel('Time(s)')
    % %     legend('show', 'Location', 'northeast')
    % end
end

%% plot the summation of real and imaginary components of channel estimations
% figure,
% subplot(211)
% plot(sumChanEstiValueReal')
% subplot(212)
% plot(sumChanEstiValueImag')

%% normalization
% repreValueFatigue = mapminmax(repreValueFatigue);

%% remove outliers in the data
% repreValueFatigue_wo_out_rev = ORM(repreValueFatigue);
% repreValueFatigue_wo_out_rev = ORM_windowed(repreValueFatigue);
repreValueFatigue_wo_out_rev = repreValueFatigue;

fatigueValues = repreValueFatigue_wo_out_rev;
% 
% figure, hold on 
% plot(repreValueFatigue)
% grid on
% title('Raw fatigue measurements')
% xlabel('Workout #')

%% plot
% figure, hold on 
% plot(fatigueValues, 'b')
grid on
title('Raw fatigue measurements')
xlabel('Workout #')
ylabel('Fatigue Level')

outlierIndex_G1 = isoutlier(fatigueValues, 'median', 'ThresholdFactor', 6);

fatigueValAfterORM_G1 = fatigueValues;

fatigueValAfterORM_G1(outlierIndex_G1) = [];

fatigueValAfterORMNormSM_G1 = smooth(fatigueValAfterORM_G1, 10)';

fatigueValAfterORMNormSM_G1(end) = [];fatigueValAfterORMNormSM_G1(1) = []; % remove the start and end point

% plot(fatigueValAfterORMNormSM_G1, 'xb')

num_G1 = 1:length(fatigueValAfterORMNormSM_G1);
[p_G1, S_G1] = polyfit(num_G1, fatigueValAfterORMNormSM_G1, 1);
[fatigueValAfterORMNormSMFit_G1, ~] = polyval(p_G1, num_G1, S_G1);
% num_G2 = 1:length(fatigueValAfterORMNormSM_G1);
% [p_G2, S_G2] = polyfit(num_G2, fatigueValAfterORMNormSM_G1, 2);
% [fatigueValAfterORMNormSMFit_G2, ~] = polyval(p_G2, num_G2, S_G2);
% num_G3 = 1:length(fatigueValAfterORMNormSM_G1);
% [p_G3, S_G3] = polyfit(num_G3, fatigueValAfterORMNormSM_G1, 3);
% [fatigueValAfterORMNormSMFit_G3, ~] = polyval(p_G3, num_G3, S_G3);

plot(fatigueValAfterORMNormSMFit_G1, 'b', 'LineWidth', 2)
% plot(fatigueValAfterORMNormSMFit_G2, 'r', 'LineWidth', 2)
% plot(fatigueValAfterORMNormSMFit_G3, 'm', 'LineWidth', 2)
grid on
title('Fitting results')
xlabel('Workout #')
ylabel('Fatigue Level')
% legend('Raw', '1st order', '2nd order', '3rd order')

end
%% normalization
% repreValueFatigue_wo_out_rev = mapminmax(repreValueFatigue_wo_out_rev);

% var(repreValueFatigue)
% figure, hold on 
% plot(repreValueFatigue)
% grid on
% title('Raw fatigue measurements')
% xlabel('Workout #')
% %
% figure, hold on
% plot(repreValueFatigue_wo_out_rev,'x' , 'LineWidth', 2.0)
% repreValueFatigue_wo_out_sm = smooth(repreValueFatigue_wo_out_rev)';
% plot(repreValueFatigue_wo_out_sm,'LineWidth', 2.0)
% num = 1:length(repreValueFatigue_wo_out_sm);
% [p, S] = polyfit(num, repreValueFatigue_wo_out_sm, 3);
% [repreValueFatigue_fit, delta] = polyval(p, num, S);
% plot(repreValueFatigue_fit, 'LineWidth', 2.0)
% legend('Raw','5-point average','curve fitting'),grid on
% % axis([0 inf 0 10]);
% divider_indices = strfind(filepath,'\');
% identifier = filepath(divider_indices(end)+1:end);
% title(['Fatigue Indicator @', identifier, ' (error = ', num2str(sum(abs(delta))) ,')']), 
% xlabel('Workout #')
% ylabel('Fatigue Level')
% % compute the slope
% slope = (repreValueFatigue_fit(end) - repreValueFatigue_fit(1))/(length(repreValueFatigue_fit) - 1);
% annotation('textbox', [.2 .2 .7 .7], 'String', ['Slope = ', num2str(slope)], 'FitBoxToText', 'on')
