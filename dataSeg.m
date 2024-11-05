% This is the data segmentation for spectrum
% parameters
%%seqTrans = [0,0,1,0,0,1,0,1,1,1,0,0,0,0,1,0,0,0,1,0,0,1,0,1,1,1];
%% 1. Process data
fs = 48000;
T = 0.001;
n = fs*T;
sequenceLength = 26;
% Get original data after demodulation
variable_name = 'raw7_0308_rr;'; % Modify the variable name here! (= =)(= =)(= =)
eval(['data_dem', '=', variable_name]);
identifier_field = '(Ruirong)';   % Modify the ID field here!
% saved_num = 120;                      % ending index indicator          need to change for each data file
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
diff_step = 5;
pca_com = score(:,component); %figure, plot(pca_com)
pca_com_diff = [];
diff_index = [];
for i = 1:diff_step:length(pca_com) - 1
    pca_com_diff = [pca_com_diff, (pca_com(i + 1) - pca_com(i))/diff_step];
    diff_index = [diff_index, i];
end
% pca_com = [pca_com(1:161); pca_com(321:391); pca_com(162:end)];
% pca_com_diff = [pca_com_diff(1:161), pca_com_diff(321:391), pca_com_diff(162:end)];
% figure, plot(pca_com), title('Channel Estimation')
% figure, plot(pca_com_diff), title('Channel Estimation Slope')
figure, plot(diff_index, pca_com_diff)

%% segment manually (comments record the end index)                         % index for step 5
seg1_diff = pca_com_diff(1:44);         % 305 227 192 200 220 178 219       61 45 38 40 44 36
seg2_diff = pca_com_diff(45:75);       % 600 437 354 370 425 338 373       120 87 71 74 85 68
seg3_diff = pca_com_diff(76:107);       % 881 608 535 524 592 506 537      176 122 107 105 118 101
seg4_diff = pca_com_diff(108:139);      % 1074 824 727 711 786 657 696      215 165 145 142 157 131
seg5_diff = pca_com_diff(140:end);      % 
% % 
% seg1 = pca_com(1:219);       
% seg2 = pca_com(220:373);       
% seg3 = pca_com(374:537);       
% seg4 = pca_com(538:696);     
% seg5 = pca_com(697:end);      

% segment equally

% seg1_diff = pca_com_diff(50: 50+floor(length(pca_com_diff)/5));
% seg2_diff = pca_com_diff(50+floor(length(pca_com_diff)/5) + 1: 51+2*floor(length(pca_com_diff)/5));
% seg3_diff = pca_com_diff(51+2*floor(length(pca_com_diff)/5) + 1: 52+3*floor(length(pca_com_diff)/5));
% seg4_diff = pca_com_diff(52+3*floor(length(pca_com_diff)/5) + 1: 53+4*floor(length(pca_com_diff)/5));
% seg5_diff = pca_com_diff(53+4*floor(length(pca_com_diff)/5) + 1: end);

% seg1 = pca_com(1: 1+floor(length(pca_com)/5));
% seg2 = pca_com(1+floor(length(pca_com)/5) + 1: 2+2*floor(length(pca_com)/5));
% seg3 = pca_com(2+2*floor(length(pca_com)/5) + 1: 3+3*floor(length(pca_com)/5));
% seg4 = pca_com(3+3*floor(length(pca_com)/5) + 1: 4+4*floor(length(pca_com)/5));
% seg5 = pca_com(4+4*floor(length(pca_com)/5) + 1: end);
% 
% seg1_diff = pca_com_diff(1: 1+floor(length(pca_com_diff)/5));
% seg2_diff = pca_com_diff(1+floor(length(pca_com_diff)/5) + 1: 2+2*floor(length(pca_com_diff)/5));
% seg3_diff = pca_com_diff(2+2*floor(length(pca_com_diff)/5) + 1: 3+3*floor(length(pca_com_diff)/5));
% seg4_diff = pca_com_diff(3+3*floor(length(pca_com_diff)/5) + 1: 4+4*floor(length(pca_com_diff)/5));
% seg5_diff = pca_com_diff(4+4*floor(length(pca_com_diff)/5) + 1: end);


% var_seg = []; var_diff_seg = [];
% mean_seg = []; mean_diff_seg = [];
% energy_seg = [];
% power_seg = [];
% ZCR = [];
% var_seg = [var_seg, var(seg1), var(seg2), var(seg3), var(seg4), var(seg5)];
var_diff_seg = [var_diff_seg, var(seg1_diff)/length(seg1_diff), var(seg2_diff)/length(seg2_diff), var(seg3_diff)/length(seg3_diff), var(seg4_diff)/length(seg4_diff), var(seg5_diff/length(seg5_diff))];
% figure, plot(var_seg), hold on
% plot(smooth(var_seg))
% title('Channel Estimation Slope Variance (Experiment 2)')
% legend('Variance', 'Variance after smoothing')
% xlabel('workout #')
% mean_seg = [mean_seg, mean(seg1), mean(seg2), mean(seg3), mean(seg4), mean(seg5)];
% mean_diff_seg = [mean_diff_seg, mean(seg1_diff), mean(seg2_diff), mean(seg3_diff), mean(seg4_diff), mean(seg5_diff)];
% F = fft(seg1);
% energy_seg = [energy_seg, sum(F.*conj(F))];
% 
% F = fft(seg2);
% energy_seg = [energy_seg, sum(F.*conj(F))];
% 
% F = fft(seg3);
% energy_seg = [energy_seg, sum(F.*conj(F))];
% 
% F = fft(seg4);
% energy_seg = [energy_seg, sum(F.*conj(F))];
% 
% F = fft(seg5);
% energy_seg = [energy_seg, sum(F.*conj(F))];
% power_seg = [power_seg, rms(seg1)^2, rms(seg2)^2, rms(seg3)^2, rms(seg4)^2, rms(seg5)^2];
% ZCR = [ZCR, mean(abs(diff(sign(seg1)))), mean(abs(diff(sign(seg2)))), mean(abs(diff(sign(seg3)))), ...
%     mean(abs(diff(sign(seg4)))), mean(abs(diff(sign(seg5))))];

%% 3. Plot figure
figure
subplot(221)
% figure, 
plot(mean_seg);
title(['Channel Estimation Mean ', identifier_field])
hold on
plot(smooth(mean_seg))
legend('Mean', 'Mean after smoothing')
xlabel('workout #')
grid on

subplot(222)
% figure, 
plot(var_seg);
title(['Channel Estimation Variance ', identifier_field])
hold on
plot(smooth(var_seg))
legend('Variance', 'Variance after smoothing')
xlabel('workout #')
grid on

subplot(223)
% figure, 
plot(mean_diff_seg);
title(['Channel Estimation Slope Mean ', identifier_field])
hold on
plot(smooth(mean_diff_seg))
legend('Slope Variance', 'Slope Mean after smoothing')
xlabel('workout #')
grid on

subplot(224)
% figure, 
plot(var_diff_seg);
title(['Channel Estimation Slope Variance ', identifier_field])
hold on
plot(smooth(var_diff_seg))
legend('Slope Variance', 'Slope Variance after smoothing')
xlabel('workout #')
grid on

%
% ------------------------------------------
%
% diff_thres = 0.45;       % Modify the threshold value here! (= =)(= =)(= =)
% pos_bigvar = [];
% i = 10;      % must have a initial offset here
% guard_period = floor(length(pca_com_diff)/10);   % variance start detection guard period value
% while i <= length(pca_com_diff)
%     if abs(pca_com_diff(i)) >= diff_thres             
%         pos_bigvar = [pos_bigvar, i]
%         i = i + guard_period;    
%     else
%         i = i + 1;
%     end
% end
% % segment the data with start detection results
% p = flipud(p);
% [p_r1,p_c1] = size(p);
% % map the difference to the spectrogram
% rat_spec_diff = p_c1/length(pca_com_diff);
% size_seg = ceil(p_c1 / 10);         % size of each segment could be calculated here
% rat_prior_post = 0.3;               % post = 1 - prior
% for i = 1:length(pos_bigvar)
%     % segment the spectrogram
%     begin_point_spec = floor(pos_bigvar(i) * rat_spec_diff) - ceil(rat_prior_post * size_seg);
%     end_point_spec = floor(pos_bigvar(i) * rat_spec_diff) + ceil((1 - rat_prior_post) * size_seg);
%     data_for_save = p(:,begin_point_spec:end_point_spec);
%     % normalize the data
%     data_norm = mapminmax(data_for_save, 0, 1);
%     figure, imagesc(db(data_norm));
%     filename = ['MUS_DATASET\', identifier_field, '_', variable_name(end-5:end-1), '\', identifier_field, '_', variable_name(end-5:end-1), '_', num2str(i + saved_num), '.csv'];
%     csvwrite(filename, data_norm);
% end
