% This is to process the acceleration data files from muscle accelerometer


filepath = uigetdir('F:\Matlab_Files\muscle_v3\accData','Select Folder for accelerometer files');
if filepath == 0
    return
end
listing = dir([filepath,'\*.csv']);
accReadings = zeros(500,3,length(listing));
for fileIndex = 1:length(listing)
    accReadings(:,:,fileIndex) = csvread([filepath, '\' , listing(fileIndex).name]);
end

%% Process
accCombine = zeros(500, length(listing));
for i = 1:length(listing)
    for j = 1:500
        accCombine(j, i) = sqrt(accReadings(j, 1, i).^2 + accReadings(j, 2, i).^2 + accReadings(j, 3, i).^2);
    end
end

% for i = 1:length(listing)
%     figure, plot(accCombine(:, i))
% end
% 
% specWoAcc = zeros(3, 300, 22);
% for i = 1:22
%     for j = 1:3
%         acousticSpec = fft(acousticData(j, :, i));
%         accSpec = fft(accCombine(:, i));
%         acousticSpecLow = acousticSpec(1:61)';
%         specDiff = acousticSpecLow - accSpec;
%         acousticSpec(1:61) = specDiff';
%         specWoAcc(j, :, i) = acousticSpec;
%     end
% end
% 
% acousticDataWoAcc = zeros(300, 22);
% temp = zeros(3, 300);
% for i = 1:22
%     for j = 1:3
%         temp(j, :) = real(ifft(specWoAcc(j, :, i)));
%     end
%     acousticDataWoAcc(:, i) = mean(temp);
% end
% figure, plot(mean(acousticDataWoAcc))
% 
% 
% colorMap = [0 0 0;
%             0 0 1;
%             0 1 0;
%             0.5 0 0.5;
%             1 0 0;
%             1 0 1;
%             0.5 0.5 0;
%             0 0 0.5;
%             0 0.5 0;
%             0 0.5 0.5;
%             0.5 0 0];

for i = 1:size(accCombine, 2)
%     acousticDataSmooth = smooth(acousticData(1, :, i));


   figure, 
%    subplot(211)
%    hold on
% %    for j = 1:6
% %        plot(acousticData(1,((j-1)*50 + 1) : j *50 ,i), 'Color', colorMap(j,:), 'LineWidth', 2)
% %    end
%    plot(acousticData(1,:,i))
%    grid on
   suptitle(['Workout ', num2str(i)])
%    xlabel('Sample Index #')
%    ylabel('Acoustic Metrics')
   
   
%    xlabel('I component')
%    ylabel('Q component')
%    legend('1-50','51-100','101-150','151-200','201-250','251-300')
    subplot(211)
    plot(accCombine(:, i))
    grid on
    xlabel('Sample Index #')
    ylabel('Accelerometer Metrics')
    
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    x0=100;
    y0=100;
    width=550;
    height=400;
    set(gcf,'position',[x0,y0,width,height])
   
%% psd
%    subplot(222)
%    [psd, f] = periodogram(acousticData(1,:,i), [], 128, 100);
% %     plot(f, 10*log10(psd));
%     plot(f(2:end), psd(2:end));
    
    subplot(212)
   [psd, f] = periodogram(accCombine(:, i), [], 128, 100);
%     plot(f, 10*log10(psd));
    plot(f(2:end), psd(2:end));
    grid on
    xlabel('Frequency')
    ylabel('PSD')
    
    
   ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    x0=100;
    y0=100;
    width=550;
    height=400;
    set(gcf,'position',[x0,y0,width,height])
end