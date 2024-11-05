% This is the data segmentation for k means clustering
% THIS SCRIPT MUST BE RUN WITH BREAKPOINTS!!!!!!!!!!!!!!!!!
% parameters
%%seqTrans = [0,0,1,0,0,1,0,1,1,1,0,0,0,0,1,0,0,0,1,0,0,1,0,1,1,1];
% k_means_tmp = [];     % THIS IS TO INITIALIZE THE K-MEANS MATRIX
function k_means_tmp = kmeans_segformat(variable, k_means_tmp)
fs = 48000;
T = 0.001;
n = fs*T;
P = 16;
L = 10;
M = [];
seqTrans = [-1,-1,1,-1,-1,1,-1,1,1,1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,-1,1,-1,1,1,1]; % for bpsk
%seqTrans = [-1,-1,1,-1,1,1,-1,1,1,1,-1,1,1,1,1,-1,-1,-1,1,-1,1,1,-1,1,1,1];
for i = 1 : L+1
    M = [seqTrans(i : i + P - 1)', M];
end
sequenceLength = 26;
% Get original data after demodulation
% variable_name = 'bpsk_0225_3_14d0p;'; % Modify the variable name here! (= =)(= =)(= =)
% eval(['data_dem', '=', variable_name, ';']);
% identifier_field = 'Xingzhe';   % Modify the ID field here!
% saved_num = 0;                      % ending index indicator          need to change for each data file
data_dem = variable;
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
%%p = pspectrum(score(:,component),38,'spectrogram', 'TimeResolution',1);
%%figure, plot((1:l_c1) * sequenceLength * T, score(:,component))
% difference
diff_step = 1;
cutoff = 50;
pca_com = score(:,component);
figure, plot(pca_com)
pca_com_diff = [];
for i = 1:diff_step:length(pca_com) - 1
    pca_com_diff = [pca_com_diff, (pca_com(i + 1) - pca_com(i))/diff_step];
end
pca_com = pca_com(cutoff:end);
pca_com = mapminmax(pca_com, 0, 1);
% pca_com = (pca_com - mean(pca_com))./std(pca_com);
diff_thres = 0.9;       % Modify the threshold value here! (= =)(= =)(= =)
diff_thres_ = 0.22;
pos_bigvar = [];
i = 1;      % must have a initial offset here
% guard_period = floor(length(pca_com_diff)/20);   % variance start detection guard period value
% forward and backward search
while i <= length(pca_com)
    if abs(pca_com(i)) >= diff_thres && (abs(pca_com_diff(i+cutoff)) >= diff_thres_|| abs(pca_com_diff(i+2+cutoff)) >= diff_thres_|| abs(pca_com_diff(i-2+cutoff)) >= diff_thres_|| abs(pca_com_diff(i+4+cutoff)) >= diff_thres_)             
        pos_bigvar = [pos_bigvar, i]
        break;
%         i = i + guard_period;    
    else
        i = i + 1;
    end
end
diff_thres = 0.9;       % Modify the threshold value here! (= =)(= =)(= =)
i = length(pca_com);
while i >= 1
    if abs(pca_com(i)) >= diff_thres && i - pos_bigvar(1) <= 80 && (abs(pca_com_diff(i+cutoff)) >= diff_thres_ || abs(pca_com_diff(i-2+cutoff)) >= diff_thres_ || abs(pca_com_diff(i+2+cutoff)) >= diff_thres_ || abs(pca_com_diff(i-4+cutoff)) >= diff_thres_ ) % 48000/26/48, 1 workout in 1 sec           
        pos_bigvar = [pos_bigvar, i]
        break;
%         i = i + guard_period;    
    else
        i = i - 1;
    end
end
% segment the data with start detection results
size_seg = ceil(length(pca_com) / 5);         % size of each segment could be calculated here
rat_prior_post = 0.2;               % post = 1 - prior
rat_prior_post_ = 0.4;               % post = 1 - prior
scaled_length = 128;               % segment length after length
for i = 1:1
    % segment the data
    begin_point_spec = pos_bigvar(i) - rat_prior_post .* size_seg;
    end_point_spec = pos_bigvar(i + 1) + rat_prior_post_ .* size_seg;
    data_for_save = pca_com(floor(begin_point_spec):ceil(end_point_spec));
%     figure, plot(data_for_save);
    % normalize the data
%     data_norm = mapminmax(data_for_save, 0, 1);
    data_norm = (data_for_save - mean(data_for_save))./std(data_for_save);
    % scale the time series to 1024
    x_before_interp = 1:length(data_norm);
    x_after_interp = linspace(1,length(data_norm),scaled_length);
    y_interp = interp1(x_before_interp,data_norm,x_after_interp);
% %     figure, plot(y,'r');
    y = y_interp;
    k_means_tmp = [k_means_tmp, y'];    % comment here to see the result first
    figure, plot(pca_com)
    hold on, plot([pos_bigvar(1), pos_bigvar(1)],[-1, 1],'r');
    plot([pos_bigvar(2), pos_bigvar(2)],[-1, 1],'r');
    figure, plot(pca_com_diff)
    figure, plot(y);
%     y_interp = reshape(y_interp, [32, 32]);
%     figure, imagesc(db(y_interp));
    % save as csv file
%     filename = ['MUS_DATASET_TIMESERIES1\', identifier_field, '_', variable_name(end-5:end-1), '\', identifier_field, '_', variable_name(end-5:end-1), '_', num2str(i + saved_num), '.csv'];
%     csvwrite(filename, y_interp);
end
% legend('1','2','3','4','5')
% figure, plot(pca_com_diff)
