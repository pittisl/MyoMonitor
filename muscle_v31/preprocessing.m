% This is to process the received data by computing its channel estimation
% This has to use the data collected by MyoTester App
function chanEstimations = preprocessing(pcm_data, maxChanTap)
    %% Parameters
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
    guardSamples = 10000;
    perTapDist = 100 * c/(2 * fc/codeUpsampleRate); % in centimeter
    load('M_mat.mat'); % load matrix M 
    load('ref_code.mat'); % load reference code file in variable 
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

    %% Compute autocorrelation to detect the start point
    corr = zeros(1, guardSamples);
    for i = 1:guardSamples
        R = corrcoef(pcm_data(i:i+length(y)-1),y);
        corr(i) = R(2,1);
    end

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
    
    %% Compute the difference h^t - h^t-1
    h_d_ = zeros(l_r1,l_c1-1);
    for i = 1 : l_c1-1
        h_d_(:,i) = h_(:,i + 1) - h_(:,i);
    end
    
    centerIcom = zeros(L+1, 1);
    centerQcom = zeros(L+1, 1);
    sumDistance = zeros(L+1, size(h_, 2));
    for i = 1 : L+1
        centerIcom(i) = mean(real(h_(i, :)));
        centerQcom(i) = mean(imag(h_(i, :)));
        [~,~,sumDistance(i, :)] = sumDistanceFunc(h_(i, :), centerIcom(i), centerQcom(i));
    end
    chanEstimations = abs(sumDistance(2:maxChanTap, :));
%     figure, plot(chanEstimations')
%     chanEstimations = abs(h_(2:maxChanTap, :));
%     chanEstimations = h_(2:maxChanTap, :);
    
    %% Plot
%     figure
%     for i = 1:maxChanTap
%         subplot(maxChanTap,1,i), hold on
%         plot(chanEstimations(i, :), 'LineWidth', 1)
%     end
end