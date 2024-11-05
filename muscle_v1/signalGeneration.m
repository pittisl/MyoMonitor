% This is to generate the signal for playing
clear;clc;

% Parameters setting
seqSymbols = [-1,-1, 1,-1,-1, ...
               1,-1, 1, 1, 1, ...
              -1,-1,-1,-1, 1, ...
              -1,-1,-1, 1,-1, ...
              -1, 1,-1, 1, 1, 1];
gapSymbols = repelem(0,24);
basebandSymbols = [seqSymbols, gapSymbols];

fs = 48000;         % sample rate
f_lo = 18000;       % lower band
f_hi = 22000;       % upper band
fc = 20000;         % central frequency
T = 20;             % sound duration

% cyclesPerCode = fs/(f_hi - f_lo);   % # of cycles every bit code 
cyclesPerCode = 4;     % IF MODIFY THIS, ref_code.mat SHOULD ALSO BE RENEWED

% Generate the symbols with fs and fc
t = 0:1/fs:1/fc*cyclesPerCode*length(basebandSymbols);   % symbols at each interval 1/fs
basebandSymbols_fs = zeros(1,round(fs/fc*cyclesPerCode*length(basebandSymbols)));
n_time = 0;
for i = 1:length(basebandSymbols)
    while n_time * 1/fs < i * 1/fc*cyclesPerCode
        n_time = n_time + 1;
        basebandSymbols_fs(n_time) = basebandSymbols(i);
    end
end
x = sqrt(2)*cos(2*pi*fc*t); % Carrier wave
y = x(1:end-1).*basebandSymbols_fs; % modulation
% figure, plot(y)

% extend to specified time duration
copyCnt = round(fs*T/length(y));
y = repmat(y, 1, copyCnt);
% figure, plot(y)

% Remove discontinuity
y_freq = fft(y);
y_freq(1:length(y_freq)*f_lo/fs) = 0;
y_freq(length(y_freq)*(1 - f_lo/fs): end) = 0; 
y_freq(round(length(y_freq)*(f_hi/fs)): round(length(y_freq)*(1-f_hi/fs))) = 0; % BPF keeping only 18-22k
y_freq = [y_freq(1:length(y_freq)/2), repelem(0,length(y) * (cyclesPerCode - 1)),y_freq((length(y_freq)/2) + 1:end)]; % padding zeros
y = real(ifft(y_freq));
% figure, plot(y)

quant=max(y)/(2^15-1);
z=round(y/quant);
%plot(z)
outint = int16(z);
%signe=uint16((sign(z)'+1)/2);
%out=[signe dec2bin(abs(z),15)];

audiowrite('sequence_bpsk_4.wav', outint, fs * cyclesPerCode);