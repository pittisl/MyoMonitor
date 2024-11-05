% This is to process the received data
clear;clc;

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
dura_stat = 3; dura_power = 3; dura_lift = 3; dura_hold = 3; dura_put = 4; % standard test timing
dura_total = dura_stat + dura_power + dura_lift + dura_hold + dura_put;
        
%% Read data from PCM file
[filename, pathname] = uigetfile('Data/*.pcm', 'Please select');
pcm_data = ReadAudioFile([pathname, filename]);
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
% figure, hold on
% for i = 1 : L+1
%     subplot(4,3,i), hold on
%         plot(h_(i, :), ...
%         'Color', colorMap(i,:), 'LineWidth', 2)
%     title([num2str(perTapDist * (i-1)), ' to ',num2str(perTapDist * i), 'cm'])
%     legend('show', 'Location', 'southeast')
% end
% suptitle('IQ traces from different channel tap')

%% Plot the real and imag component
% h_abs = zeros(l_r1, l_c1);
% figure, hold on
% for i = 1 : L+1
%     subplot(4,3,i), plot(linspace(0, l_c1 * frameLength * samplesPerCode / fs, l_c1)+startPoint/fs, ...
%         (real(h_(i, :))/(2*pi))*(c/fc)*0.5*100, ...
%         'Color', colorMap(i,:), 'DisplayName', [num2str(perTapDist * (i-1)), ' to ',num2str(perTapDist * i), 'cm'], ...
%         'LineWidth', 2)
% %     h_abs(i, :) = (smooth(phase(h_(i, :)))/(2*pi))*(c/fc)*0.5*100;
%     xlabel('Time(s)')
%     legend('show', 'Location', 'southwest')
% %     axis([-inf inf -1 1]);
% end
% suptitle('Real component of channel estimation')
% 
% figure, hold on
% for i = 1 : L+1
%     subplot(4,3,i), plot(linspace(0, l_c1 * frameLength * samplesPerCode / fs, l_c1)+startPoint/fs, ...
%         (imag(h_(i, :))/(2*pi))*(c/fc)*0.5*100, ...
%         'Color', colorMap(i,:), 'DisplayName', [num2str(perTapDist * (i-1)), ' to ',num2str(perTapDist * i), 'cm'], ...
%         'LineWidth', 2)
%     xlabel('Time(s)')
%     legend('show', 'Location', 'southwest')
%     axis([-inf inf -1 1]);
% end
% suptitle('Imag component of channel estimation')

%% FFT on channel estimation
h_spec = fft(h_');  % the transpose'
% figure, hold on
% for i = 1 : L+1
%     subplot(4,3,i), plot(linspace(1, fc/codeUpsampleRate/frameLength, l_c1-1), ...
%         abs(h_spec(2:end, i)), ...
%         'Color', colorMap(i,:), 'DisplayName', [num2str(perTapDist * (i-1)), ' to ',num2str(perTapDist * i), 'cm'], ...
%         'LineWidth', 1)
%     xlabel('Freq(Hz)')
%     legend('show', 'Location', 'northwest')
%     axis([-inf inf 0 500]);
% end
% suptitle('Spectrum Amplitude of Channel Estimation')

%% Take out spectrum of 8-12 Hz
% sampleRate = fc/codeUpsampleRate/frameLength;
% spec_start = l_c1 * (19/sampleRate);
% spec_end = l_c1 * (28/sampleRate);
% % filter out other frequency
% h_spec(1:floor(spec_start - 1), :) = 0;
% h_spec(ceil(spec_end + 1): end, :) = 0;
% h_filtered = ifft(h_spec);

%% plot the complex value of h_
figure, hold on
for i = 1 : L+1
    subplot(4,3,i), hold on
        plot(smooth(abs(h_(i, :))), ...
        'Color', colorMap(i,:), 'LineWidth', 2)
    title([num2str(perTapDist * (i-1)), ' to ',num2str(perTapDist * i), 'cm'])
%     legend('show', 'Location', 'southeast')
end
suptitle('IQ traces from different channel tap')


% % h_filtered_abs = abs(h_filtered);
% % sumChan1to3 = sum(h_filtered_abs(:,1:3),2);% sum beginning 3 channels
% % h_filtered_abs_sm = smooth(sumChan1to3, 50);
% % % figure,plot(h_filtered_abs_sm);
% % h_lift = h_filtered_abs_sm(floor(length(h_filtered_abs_sm)*6/16):ceil(length(h_filtered_abs_sm)*9/16));
% % figure,plot(h_lift)
% % sum(h_lift)

% figure, hold on
% for i = 1 : L+1
%     subplot(4,3,i), plot(linspace(8, 12, ceil(spec_end) - floor(spec_start) + 1), ...
%         smooth(abs(h_spec(floor(spec_start):ceil(spec_end), i))), ...
%         'Color', colorMap(i,:), ...
%         'LineWidth', 1)
%     title([num2str(perTapDist * (i-1)), ' to ',num2str(perTapDist * i), 'cm']);
%     xlabel('Freq(Hz)')
% %     legend('show', 'Location', 'northwest')
%     axis([-inf inf 0 100]);
% end
% suptitle('Spectrum Amplitude of Channel Estimation (8-12Hz)')

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
centerIcom = zeros(L+1, 4);
centerQcom = zeros(L+1, 4);
sumDistance = zeros(L+1, 4);
% figure,
for i = 1 : L+1
%     subplot(4,3,i), hold on
%     plot(h_(i, period_stat), ...
%         'LineWidth', 1, 'DisplayName', 'static')
%     plot(h_(i, period_power), ...
%         'LineWidth', 1, 'DisplayName', 'generate')
%     plot(h_(i, period_lift'), ...
%         'LineWidth', 1, 'DisplayName', 'lift-up')
%     plot(h_(i, period_hold), ...
%         'LineWidth', 1, 'DisplayName', 'hold')
% %     plot(h_(i, period_put), ...
% %         'LineWidth', 1, 'DisplayName', 'lift-down')
%     title([num2str(perTapDist * (i-1)), ' to ',num2str(perTapDist * i), 'cm'])
%     legend('show', 'Location', 'northeast')
    % calculate the centroid of the point set in each workout section and
    % sum the distance of all the points
    centerIcom(i,:) = [mean(real(h_(i, period_stat))), mean(real(h_(i, period_power))), ...
        mean(real(h_(i, period_lift))), mean(real(h_(i, period_hold)))];
    centerQcom(i,:) = [mean(imag(h_(i, period_stat))), mean(imag(h_(i, period_power))), ...
        mean(imag(h_(i, period_lift))), mean(imag(h_(i, period_hold)))];
    sumDistance(i,:) = [sumDistanceFunc(h_(i, period_stat), centerIcom(i,1), centerQcom(i,1)), ...
        sumDistanceFunc(h_(i, period_power), centerIcom(i,2), centerQcom(i,2)), ...
        sumDistanceFunc(h_(i, period_lift), centerIcom(i,3), centerQcom(i,3)), ...
        sumDistanceFunc(h_(i, period_hold), centerIcom(i,4), centerQcom(i,4))];
end
% suptitle('IQ traces from different channel tap with workout label')

% figure, 
% for i = 1 : L+1
%     subplot(4,3,i), hold on
%     plot(h_filtered(period_stat, i), ...
%         'LineWidth', 1, 'DisplayName', 'static')
% %     plot(h_filtered(period_power, i), ...
% %         'LineWidth', 1, 'DisplayName', 'generate')
%     plot(h_filtered(period_lift, i'), ...
%         'LineWidth', 1, 'DisplayName', 'lift-up')
% %     plot(h_filtered(period_hold, i), ...
% %         'LineWidth', 1, 'DisplayName', 'hold')
% %     plot(h_(i, period_put), ...
% %         'LineWidth', 1, 'DisplayName', 'lift-down')
%     title([num2str(perTapDist * (i-1)), ' to ',num2str(perTapDist * i), 'cm'])
%     legend('show', 'Location', 'northeast')
%     axis([-1 1 -1 1]);
%     % calculate the centroid of the point set in each workout section and
%     % sum the distance of all the points
% end
% suptitle('IQ traces from different channel tap with workout label')

%% plot the degree of concentration
% figure,hold on
% for i = 1 : L+1
%     plot(sumDistance(i,:), 'Color', colorMap(i,:), ...
%         'DisplayName', [num2str(perTapDist * (i-1)), ' to ',num2str(perTapDist * i), 'cm'], ...
%         'LineWidth', 2)
%     legend('show', 'Location', 'northeast')
% end
% title('IQ traces level of concentration at each workout stage')

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
% figure, hold on
% for i = 1 : L+1
%     subplot(4,3,i), plot(linspace(0, l_c1 * frameLength * samplesPerCode / fs, l_c1)+startPoint/fs, ...
%         smooth(angle(h_(i, :))), ...
%         'Color', colorMap(i,:), 'DisplayName', [num2str(perTapDist * (i-1)), ' to ',num2str(perTapDist * i), 'cm'], ...
%         'LineWidth', 2)
%     xlabel('Time(s)')
%     legend('show', 'Location', 'southwest')
% end
% suptitle('Raw Phase Channel Tap')
% xlabel('Time(s)')
% legend('show', 'Location', 'southwest')

%% Plot amplitude of channel estimation
figure, hold on
for i = 1 : L+1
    subplot(4,3,i), plot(linspace(0, l_c1 * frameLength * samplesPerCode / fs, l_c1)+startPoint/fs, ...
        smooth(abs(h_(i, :))), ...
        'Color', colorMap(i,:), 'DisplayName', [num2str(perTapDist * (i-1)), ' to ',num2str(perTapDist * i), 'cm'], ...
        'LineWidth', 2)
    xlabel('Time(s)')
    legend('show', 'Location', 'southwest')
end
suptitle('Raw Amplitude Channel Tap')
xlabel('Time(s)')
legend('show', 'Location', 'northeast')

%% Compute the difference h^t - h^t-1
h_d_ = zeros(l_r1,l_c1-1);
for i = 1 : l_c1-1
    h_d_(:,i) = h_(:,i + 1) - h_(:,i);
end

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