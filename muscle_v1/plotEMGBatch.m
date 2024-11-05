% This is to plot the emg signal IN FOLDER
% clc; clear;
function fatigueEMG = plotEMGBatch()
%% Read EMG data from text files
% [rawEMG_filename, rawEMG_pathname] = uigetfile('EMG/*.txt', 'Please select raw EMG file');
% rawEMG_data = textread([rawEMG_pathname, rawEMG_filename], '%f');
% figure, 
% subplot(211), 
% plot(linspace(0, 16, length(rawEMG_data)), rawEMG_data), 
% hold on
% % add timeline
% plot([3, 3],[min(rawEMG_data), max(rawEMG_data)],'r');
% plot([6, 6],[min(rawEMG_data), max(rawEMG_data)],'r');
% plot([9, 9],[min(rawEMG_data), max(rawEMG_data)],'g');
% plot([12, 12],[min(rawEMG_data), max(rawEMG_data)],'g');
% xlabel('Time(s)')
% title('Raw EMG signal for a single workout')
filepath = uigetdir('EMG', 'Please select averaged EMG folder');
listing = dir([filepath,'\*.txt']);
fatigueEMG = zeros(1,length(listing));
for fileIndex = 1:length(listing)
    averEMG_data = textread([filepath, '\', listing(fileIndex).name], '%f');
    transPoints = transPointsIdentifyEMG(averEMG_data);
    fatigueEMG(fileIndex) = sum(averEMG_data(transPoints(1):transPoints(4)));
end

% fatigueEMG_sm = smooth(fatigueEMG)';

% num = 1:length(fatigueEMG_sm);
% [p, S] = polyfit(num, fatigueEMG_sm, 1);
% [fatigueEMG_fit, delta] = polyval(p, num, S);
% error=sum(abs(delta));

figure
hold on
grid on
plot(fatigueEMG);
end
% plot(fatigueEMG,'x', 'LineWidth', 2.0);
% plot(fatigueEMG_sm, 'LineWidth', 2.0)
% plot(fatigueEMG_fit, 'LineWidth', 2.0);
% legend('Raw','5-point average','curve fitting')
% xlabel('Workout #')
% ylabel('Fatigue Level')
% title('Fatigue from EMG Signal')
% 
% % compute the slope
% slope = (fatigueEMG_fit(end) - fatigueEMG_fit(1))/(length(fatigueEMG_fit) - 1);
% annotation('textbox', [.2 .2 .7 .7], 'String', ['Slope = ', num2str(slope)], 'FitBoxToText', 'on')