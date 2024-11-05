% This is to plot multiple groups results for EMG signals
clear; clc;
fatigueValuesEMG_G1 = plotEMGBatch;
fatigueValuesEMG_G2 = plotEMGBatch;
fatigueValuesEMG_G3 = plotEMGBatch;

%% Raw data plot
figure, hold on 
plot(fatigueValuesEMG_G1, 'b')
plot(fatigueValuesEMG_G2, 'r')
plot(fatigueValuesEMG_G3, 'm')
grid on
title('Raw fatigue measurements from EMG')
xlabel('Workout #')
ylabel('Fatigue Level')
legend('Group 1','Group 2','Group 3')


outlierIndex_G1 = isoutlier(fatigueValuesEMG_G1, 'median');
outlierIndex_G2 = isoutlier(fatigueValuesEMG_G2, 'median');
outlierIndex_G3 = isoutlier(fatigueValuesEMG_G3, 'median');
fatigueValAfterORMEMG_G1 = fatigueValuesEMG_G1;
fatigueValAfterORMEMG_G2 = fatigueValuesEMG_G2;
fatigueValAfterORMEMG_G3 = fatigueValuesEMG_G3;
fatigueValAfterORMEMG_G1(outlierIndex_G1) = [];
fatigueValAfterORMEMG_G2(outlierIndex_G2) = [];
fatigueValAfterORMEMG_G3(outlierIndex_G3) = [];


fatigueValAfterORMSM_EMG_G1 = smooth(fatigueValAfterORMEMG_G1, 10)';
fatigueValAfterORMSM_EMG_G2 = smooth(fatigueValAfterORMEMG_G2, 10)';
fatigueValAfterORMSM_EMG_G3 = smooth(fatigueValAfterORMEMG_G3, 10)';
