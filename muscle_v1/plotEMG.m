% This is to plot the emg signal
% clc; clear;
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
[averEMG_filename, averEMG_pathname] = uigetfile('EMG/*.txt', 'Please select averaged EMG file');
averEMG_data = textread([averEMG_pathname, averEMG_filename], '%f');
figure,
% subplot(212),
plot(linspace(0, 10, length(averEMG_data)), averEMG_data)
hold on
% add adaptive timeline
transPoints = transPointsIdentifyEMG(averEMG_data);
transPoints_l = (transPoints/length(averEMG_data))*10;
plot([transPoints_l(1), transPoints_l(1)],[min(averEMG_data), max(averEMG_data)],'r');
% plot([transPoints_l(2), transPoints_l(2)],[min(averEMG_data), max(averEMG_data)],'r');
% plot([transPoints_l(3), transPoints_l(3)],[min(averEMG_data), max(averEMG_data)],'r');
plot([transPoints_l(4), transPoints_l(4)],[min(averEMG_data), max(averEMG_data)],'r');
% add fixed timeline
plot([3, 3],[min(averEMG_data), max(averEMG_data)],'g');
% plot([5, 5],[min(averEMG_data), max(averEMG_data)],'r');
% plot([8, 8],[min(averEMG_data), max(averEMG_data)],'g');
plot([10, 10],[min(averEMG_data), max(averEMG_data)],'g'); 
xlabel('Time(s)')
title('Processed EMG signal for a single workout')
sum(averEMG_data)