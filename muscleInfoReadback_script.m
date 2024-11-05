% This is the script to process the PCM file from phone

% Create reference signal
% create sequence to transmit 26-bit GSM training sequence
seqTrans = [0,0,1,0,0,1,0,1,1,1,0,0,0,0,1,0,0,0,1,0,0,1,0,1,1,1];
sequenceLength = 26;
% ASK carrier wave?
% modified with BPSK
fs = 48000;
T = 0.001;
n = fs*T;
f = 18000;
y = sin(2*pi*f*T*linspace(0,1,n));
% % tuper window
% t5_w = tukeywin(n,0.8);
% y = y.*t5_w';
% % tuper window end
dura = 60;
% figure, plot(y)
% convert sequence to raw modulated signal
yy = [];
for i = 1 : length(seqTrans)
    if seqTrans(i) == 1
        yy = [yy, y];
    else
%         yy = [yy, zeros(1, fs*T)];
        yy = [yy, -y];
    end
end
% figure, plot(yy)
raw_wav = zeros(1, fs * dura);
for i = 1 : fs * dura
    if mod(i,length(yy)) == 0
        raw_wav(i) = yy(1248);
    else
        raw_wav(i) = yy(mod(i, length(yy)));
    end
end

% read PCM file
pcmfile = uigetfile('.pcm','Please select PCM');
pcm_data = ReadAudioFile(pcmfile); % ReadAudioFile should be under same directory
% identify the start point
corr = [];
for i = 1:100000
    R = corrcoef(pcm_data(i:i+1247),yy);
    corr(i) = R(2,1);
end
for i = 1:length(corr)
    if corr(i) >= 0.75
        startPoint = i;
        break;
    end
end

% demodulate
pcm_data_cor = pcm_data(startPoint:end);
y_c = raw_wav(1:length(pcm_data_cor));
y_mul = y_c .* pcm_data_cor';
% % ---------------
% % y_mul_spec = fft(y_mul);
% % y_mul_spec(50:end - 50) = 0;
% % y_dem = ifft(y_mul_spec);
% % figure,plot(real(y_dem))
% % ---------------
% % figure, plot(y_mul)
y_filter = fir1(16, 1000*2/48000);  %LPF
y_dem = filter(y_filter, 1, y_mul);
% figure,plot(y_dem)
% figure,plot((1:length(y_dem))/(fs*T),y_dem)
y_dem_pad = [zeros(1, startPoint - 1), y_dem];
% figure,plot((1:(length(y_dem_pad)))/fs,y_dem_pad/max(y_dem_pad)*2)

% % --- generate reference signal ---
% yr = [];
% for i = 1 : length(seqTrans)
%     if seqTrans(i) == 1
%         yr = [yr, ones(1, fs*T)];
%     else
%         yr = [yr, zeros(1, fs*T)];
%     end
% end
% % figure, plot(yr)
% yr_all = zeros(1, length(y_dem));
% for i = 1 : length(y_dem)
%     if mod(i,length(yr)) == 0
%         yr_all(i) = yr(1248);
%     else
%         yr_all(i) = yr(mod(i, length(yr)));
%     end
% end
% %figure, 
% yr_all_pad = [zeros(1, startPoint - 1), yr_all];
% hold on, plot((1:(length(y_dem_pad)))/fs, yr_all_pad)


% --- estimate the channel ---
% create M
P = 16;
L = 10;
M = [];
seqTrans = [-1,-1,1,-1,-1,1,-1,1,1,1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,-1,1,-1,1,1,1]; % for bpsk
%seqTrans = [-1,-1,1,-1,1,1,-1,1,1,1,-1,1,1,1,1,-1,-1,-1,1,-1,1,1,-1,1,1,1];
for i = 1 : L+1
    M = [seqTrans(i : i + P - 1)', M];
end

% % edit here to save the signal after demodulation
% eval([filename, '=', 'y_dem;']);

% get demodulation data
data_dem = y_dem;
% Get channel estimation vector
data_dem_scaled = data_dem/max(data_dem)*2;
data_bit_means = zeros(1, floor(length(data_dem_scaled)/n));   % bit mean value
data_bit_means_index = zeros(1, floor(length(data_dem_scaled)/n));
for i = 1 : length(data_dem_scaled)/n
    data_bit_means(i) = mean(data_dem_scaled((i-1)*n + 1:i * n));
    data_bit_means_index(i) = (i-1)*n + 1;
end
%%figure, plot(data_bit_means)
% take out the P length bit value  y = h * x + n
y_matrix = [];
for i = 1 : length(data_bit_means)/sequenceLength
    y_matrix = [y_matrix, data_bit_means((i-1)*sequenceLength + 1:(i-1)*sequenceLength + P)'];
end
[l_r,l_c] = size(y_matrix);
% y = h * x + n
h = [];
for i = 1:l_c
    h = [h, (M' * M) \ (M' * y_matrix(:,i))];
end
% normalization
[l_r1,l_c1] = size(h);
h_n = mapminmax(h);
% PCA
[coeff,score,latent] = pca(h_n');
component = 1;
%%figure, pspectrum(score(:,component),38,'spectrogram', 'TimeResolution',1);
% p = pspectrum(score(:,component),38,'spectrogram', 'TimeResolution',1);
%%figure, plot((1:l_c1) * sequenceLength * T, score(:,component))
% difference
diff_step = 1;
pca_com = score(:,component);
pca_com_diff = [];
for i = 1:diff_step:length(pca_com) - 1
    pca_com_diff = [pca_com_diff, (pca_com(i + 1) - pca_com(i))/diff_step];
end
% figure, plot(pca_com)
% figure, plot(pca_com_diff)

% % calculate mean and variance for each workout
% seg1 = pca_com(50: 50+floor(length(pca_com)/5));
% seg2 = pca_com(50+floor(length(pca_com)/5) + 1: 51+2*floor(length(pca_com)/5));
% seg3 = pca_com(51+2*floor(length(pca_com)/5) + 1: 52+3*floor(length(pca_com)/5));
% seg4 = pca_com(52+3*floor(length(pca_com)/5) + 1: 53+4*floor(length(pca_com)/5));
% seg5 = pca_com(53+4*floor(length(pca_com)/5) + 1: end);
% 
% % var_seg = [];
% % mean_seg = [];
% var_seg = [var_seg, var(seg1), var(seg2), var(seg3), var(seg4), var(seg5)];
% mean_seg = [mean_seg, mean(seg1), mean(seg2), mean(seg3), mean(seg4), mean(seg5)];

% read txt to get orientation readings
orientationFile = uigetfile('.txt', 'Please select orientation record');
orientation_readings = textread(orientationFile, '%f');
timestamp = orientation_readings(1:4:end);
orientation_1 = orientation_readings(2:4:end);
orientation_2 = orientation_readings(3:4:end);
orientation_3 = orientation_readings(4:4:end);
% read txt to get acceleration readings
accelerationFile = uigetfile('.txt', 'Please select acceleration record');
acceleration_readings = textread(accelerationFile, '%f');
timestamp_2 = acceleration_readings(1:4:end);
acceleration_1 = acceleration_readings(2:4:end);
acceleration_2 = acceleration_readings(3:4:end);
acceleration_3 = acceleration_readings(4:4:end);
% PCM data after preprocess is in pca_com_diff
% Start point is in startPoint
figure,
subplot(411)
plot(pca_com)
title('Channel Estimation')
subplot(412)
plot(pca_com_diff)
title('Slope of Channel Estimation')
subplot(413)
plot(timestamp + startPoint/48000, orientation_1), hold on
plot(timestamp, orientation_2),
plot(timestamp, orientation_3),
legend('Angle of cellphone y with north', 'Angle of tilting x', 'Angle of tilting z')
title('Phone orientation readings (degree)')
subplot(414)
plot(timestamp + startPoint/48000, acceleration_1), hold on
plot(timestamp, acceleration_2),
plot(timestamp, acceleration_3),
legend('x', 'y', 'z')
title('Phone accleration readings')