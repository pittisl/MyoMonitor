% create sequence to transmit 26-bit GSM training sequence
seqTrans = [0,0,1,0,0,1,0,1,1,1,0,0,0,0,1,0,0,0,1,0,0,1,0,1,1,1];
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

% figure,plot(raw_wav)

% --- generate ---
quant=max(raw_wav)/(2^15-1);
z=round(raw_wav/quant);
%plot(z)
outint = int16(z);
%signe=uint16((sign(z)'+1)/2);
%out=[signe dec2bin(abs(z),15)];

audiowrite('myfile_bpsk_gsm.wav', outint, fs);

% % --- read pcm file ---
filename = 'raw7_0308_rr';
pcm_data = ReadAudioFile([filename, '.pcm']);
%figure, plot(pcm_data);
% h_filter = fir1(32, 15000*2/48000, 'high');
% pcm_data = filter(h_filter, 1, pcm_data);

% % --- calculate correlation ---
corr = [];
for i = 1:100000
    R = corrcoef(pcm_data(i:i+1247),yy);
    corr(i) = R(2,1);
end
% figure,plot(corr)

for i = 1:length(corr)
    if corr(i) >= 0.5
        startPoint = i;
        break;
    end
end

% startPoint = 7977;

% % --- demodulate ---
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
figure,plot((1:(length(y_dem_pad)))/fs,y_dem_pad/max(y_dem_pad)*2)

% --- generate reference signal ---
yr = [];
for i = 1 : length(seqTrans)
    if seqTrans(i) == 1
        yr = [yr, ones(1, fs*T)];
    else
        yr = [yr, zeros(1, fs*T)];
    end
end
% figure, plot(yr)
yr_all = zeros(1, length(y_dem));
for i = 1 : length(y_dem)
    if mod(i,length(yr)) == 0
        yr_all(i) = yr(1248);
    else
        yr_all(i) = yr(mod(i, length(yr)));
    end
end
%figure, 
yr_all_pad = [zeros(1, startPoint - 1), yr_all];
hold on, plot((1:(length(y_dem_pad)))/fs, yr_all_pad)


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

% edit here to save the signal after demodulation
eval([filename, '=', 'y_dem;']);

% edit here to choose different demodulated signals (-_-)(-_-)(-_-)(-_-)(-_-)(-_-)(-_-)(-_-)
y_dem = bpsk_0118_1_4p;

% get y_matrix
y_dem_scaled = y_dem/max(y_dem)*2;
y_bit_means = zeros(1, floor(length(y_dem_scaled)/n));
y_bit_means_index = zeros(1, floor(length(y_dem_scaled)/n));
for i = 1 : length(y_dem_scaled)/n
    y_bit_means(i) = mean(y_dem_scaled((i-1)*n + 1:i * n));
    y_bit_means_index(i) = (i-1)*n + 1;
end
% figure, 
% hold on, 
% plot(y_bit_means_index, y_bit_means);

% take out the P length bit value
y_matrix = [];
for i = 1 : length(y_bit_means)/length(seqTrans)
    y_matrix = [y_matrix, y_bit_means((i-1)*length(seqTrans) + 1:(i-1)*length(seqTrans) + P)'];
end

[l_r,l_c] = size(y_matrix);
h = [];
for i = 1:l_c
    h = [h, (M' * M) \ (M' * y_matrix(:,i))];
end

[l_r1,l_c1] = size(h);
figure,plot((1:l_c1) * length(seqTrans) * T + startPoint/fs, h')
hold on, grid on
% modify the time stamp here
P1 = 3.206; P2 = 6.4; P3 = 12.2; P4 = 15.018;
plot([P1, P1],[-0.2, 0.2],'r');
plot([P2, P2],[-0.2, 0.2],'r');
plot([P3, P3],[-0.2, 0.2],'r');
plot([P4, P4],[-0.2, 0.2],'r');
legend('Entry 1', 'Entry 2', 'Entry 3', 'Entry 4', 'Entry 5', 'Entry 6', 'Entry 7', 'Entry 8', 'Entry 9', 'Entry 10', 'Entry 11')
xlabel('Time(s)')
title('Result with object 3 (heaviest) - trial 1')

figure, 
for i = 1 : L+1
    subplot(11,1,i), plot((1:l_c1) * length(seqTrans) * T + startPoint/fs, h(i , :));
    hold on,
    plot([P1, P1],[-0.2, 0.5],'r');
    plot([P2, P2],[-0.2, 0.5],'r');
    plot([P3, P3],[-0.2, 0.5],'r');
    plot([P4, P4],[-0.2, 0.5],'r');
end

% --- PCA ---
% normalize
h_n = mapminmax(h);
[coeff,score,latent] = pca(h_n');
figure, 
% subplot(121),
plot((1:l_c1) * length(seqTrans) * T + startPoint/fs, score(:,1:11))
% subplot(122)
hold on, grid on
plot((1:l_c1) * length(seqTrans) * T + startPoint/fs, score(:,1))
plot([P1, P1],[-1, 1],'r');
plot([P2, P2],[-1, 1],'r');
plot([P3, P3],[-1, 1],'r');
plot([P4, P4],[-1, 1],'r');
component = 1;
figure, pspectrum(score(:,component),38,'spectrogram', 'TimeResolution',1);
figure, plot((1:l_c1) * length(seqTrans) * T + startPoint/fs, score(:,component))
% spectrogram power average
p = pspectrum(score(:,component),38,'spectrogram', 'TimeResolution',1);
% p_mean_col = mean(p);
p_mean_col = mean(p(201:end,:));
figure, plot(p_mean_col)

% % % calculate the summation of L-1 entries
% % h_sum = zeros(1,l_c1);
% % for i = 1:l_c1
% %     h_sum(i) = sum(h(1:l_r1, i));
% % end
% % figure, plot((1:l_c1) * length(seqTrans) * T, h_sum)
