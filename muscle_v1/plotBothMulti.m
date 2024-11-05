% This is to plot multiple groups results For BOTH EMG AND ACOUSTICclear; clc;
clear; clc;
fatigueValues_G1 = signalProcessBatch;
fatigueValues_G2 = signalProcessBatch;
fatigueValues_G3 = signalProcessBatch;

fatigueValuesEMG_G1 = plotEMGBatch;
fatigueValuesEMG_G2 = plotEMGBatch;
fatigueValuesEMG_G3 = plotEMGBatch;

%% Raw data plot for Acoustic
figure,
subplot(231),hold on 
plot(fatigueValues_G1, 'b')
plot(fatigueValues_G2, 'r')
plot(fatigueValues_G3, 'm')
grid on
title('Raw fatigue measurements from acoustic')
xlabel('Workout #')
ylabel('Fatigue Level')
legend('Group 1','Group 2','Group 3')

%% Raw data plot
subplot(234), hold on 
plot(fatigueValuesEMG_G1, 'b')
plot(fatigueValuesEMG_G2, 'r')
plot(fatigueValuesEMG_G3, 'm')
grid on
title('Raw fatigue measurements from EMG')
xlabel('Workout #')
ylabel('Fatigue Level')
legend('Group 1','Group 2','Group 3')

%% outlier
outlierIndex_G1 = isoutlier(fatigueValues_G1, 'median');
outlierIndex_G2 = isoutlier(fatigueValues_G2, 'median');
outlierIndex_G3 = isoutlier(fatigueValues_G3, 'median');
fatigueValAfterORM_G1 = fatigueValues_G1;
fatigueValAfterORM_G2 = fatigueValues_G2;
fatigueValAfterORM_G3 = fatigueValues_G3;
fatigueValAfterORM_G1(outlierIndex_G1) = [];
fatigueValAfterORM_G2(outlierIndex_G2) = [];
fatigueValAfterORM_G3(outlierIndex_G3) = [];

outlierIndex_G1 = isoutlier(fatigueValuesEMG_G1, 'median');
outlierIndex_G2 = isoutlier(fatigueValuesEMG_G2, 'median');
outlierIndex_G3 = isoutlier(fatigueValuesEMG_G3, 'median');
fatigueValAfterORMEMG_G1 = fatigueValuesEMG_G1;
fatigueValAfterORMEMG_G2 = fatigueValuesEMG_G2;
fatigueValAfterORMEMG_G3 = fatigueValuesEMG_G3;
fatigueValAfterORMEMG_G1(outlierIndex_G1) = [];
fatigueValAfterORMEMG_G2(outlierIndex_G2) = [];
fatigueValAfterORMEMG_G3(outlierIndex_G3) = [];

%% smooth
fatigueValAfterORMSM_G1 = smooth(fatigueValAfterORM_G1, 10)';
fatigueValAfterORMSM_G2 = smooth(fatigueValAfterORM_G2, 10)';
fatigueValAfterORMSM_G3 = smooth(fatigueValAfterORM_G3, 10)';

fatigueValAfterORMSM_EMG_G1 = smooth(fatigueValAfterORMEMG_G1, 10)';
fatigueValAfterORMSM_EMG_G2 = smooth(fatigueValAfterORMEMG_G2, 10)';
fatigueValAfterORMSM_EMG_G3 = smooth(fatigueValAfterORMEMG_G3, 10)';

%% Nonlinear regression
% num_G1 = 1:length(fatigueValAfterORMSM_G1);
% num_G2 = 1:length(fatigueValAfterORMSM_G2);
% num_G3 = 1:length(fatigueValAfterORMSM_G3);
% 
% X_G1 = num_G1;
% Y = fatigueValAfterORMSM_G1;
% modelfun = @(b, x)(b(1).*x.^3 + b(2).*x.^2 + b(3).*x + b(4)); %+ b(4)*exp(b(5).*x));
% beta0 = [0.1 0.2 0.3 0.4 ];
% 
% % opts = statset('nlinfit');
% % opts.RobustWgtFun = 'bisquare';
% 
% beta_G1 = nlinfit(X_G1,Y,modelfun,beta0, opts);
% X_G2 = num_G2;
% Y = fatigueValAfterORMSM_G2;
% beta_G2 = nlinfit(X_G2,Y,modelfun,beta0, opts);
% X_G3 = num_G3;
% Y = fatigueValAfterORMSM_G3;
% beta_G3 = nlinfit(X_G3,Y,modelfun,beta0, opts);
% 
% YY_G1 = modelfun(beta_G1, X_G1);
% YY_G2 = modelfun(beta_G2, X_G2);
% YY_G3 = modelfun(beta_G3, X_G3);
% 
% subplot(222), hold on
% plot(YY_G1, 'b', 'LineWidth', 2)
% plot(YY_G2, 'r', 'LineWidth', 2)
% plot(YY_G3, 'm', 'LineWidth', 2)
% grid on
% 
% 
% %
% num_G1 = 1:length(fatigueValAfterORMSM_EMG_G1);
% num_G2 = 1:length(fatigueValAfterORMSM_EMG_G2);
% num_G3 = 1:length(fatigueValAfterORMSM_EMG_G3);
% 
% X_G1 = num_G1;
% Y = fatigueValAfterORMSM_EMG_G1;
% modelfun = @(b, x)(b(1).*x.^3 + b(2).*x.^2 + b(3).*x + b(4)); %+ b(4)*exp(b(5).*x));
% beta0 = [0.1 0.2 0.3 0.4 ];
% 
% % opts = statset('nlinfit');
% % opts.RobustWgtFun = 'bisquare';
% 
% betaEMG_G1 = nlinfit(X_G1,Y,modelfun,beta0, opts);
% X_G2 = num_G2;
% Y = fatigueValAfterORMSM_EMG_G2;
% betaEMG_G2 = nlinfit(X_G2,Y,modelfun,beta0, opts);
% X_G3 = num_G3;
% Y = fatigueValAfterORMSM_EMG_G3;
% betaEMG_G3 = nlinfit(X_G3,Y,modelfun,beta0, opts);
% 
% YY_G1 = modelfun(betaEMG_G1, X_G1);
% YY_G2 = modelfun(betaEMG_G2, X_G2);
% YY_G3 = modelfun(betaEMG_G3, X_G3);
% 
% subplot(224), hold on
% plot(YY_G1, 'b', 'LineWidth', 2)
% plot(YY_G2, 'r', 'LineWidth', 2)
% plot(YY_G3, 'm', 'LineWidth', 2)
% grid on

%% 2nd order Curve Fitting
num_G1 = 1:length(fatigueValAfterORMSM_G1);
[p_G1, S_G1] = polyfit(num_G1, fatigueValAfterORMSM_G1, 2);
[fatigueValAfterORMNormSMFit_G1, delta_G1] = polyval(p_G1, num_G1, S_G1);
num_G2 = 1:length(fatigueValAfterORMSM_G2);
[p_G2, S_G2] = polyfit(num_G2, fatigueValAfterORMSM_G2, 2);
[fatigueValAfterORMNormSMFit_G2, delta_G2] = polyval(p_G2, num_G2, S_G2);
num_G3 = 1:length(fatigueValAfterORMSM_G3);
[p_G3, S_G3] = polyfit(num_G3, fatigueValAfterORMSM_G3, 2);
[fatigueValAfterORMNormSMFit_G3, delta_G3] = polyval(p_G3, num_G3, S_G3);
subplot(232), hold on
plot(fatigueValAfterORMNormSMFit_G1, 'b', 'LineWidth', 2)
plot(fatigueValAfterORMNormSMFit_G2, 'r', 'LineWidth', 2)
plot(fatigueValAfterORMNormSMFit_G3, 'm', 'LineWidth', 2)
grid on
title('2nd order fitting on acoustic results')
xlabel('Workout #')
ylabel('Fatigue Level')


num_G1 = 1:length(fatigueValAfterORMSM_EMG_G1);
[p_G1, S_G1] = polyfit(num_G1, fatigueValAfterORMSM_EMG_G1, 2);
[fatigueValAfterORMNormSMFit_G1, delta_G1] = polyval(p_G1, num_G1, S_G1);
num_G2 = 1:length(fatigueValAfterORMSM_EMG_G2);
[p_G2, S_G2] = polyfit(num_G2, fatigueValAfterORMSM_EMG_G2, 2);
[fatigueValAfterORMNormSMFit_G2, delta_G2] = polyval(p_G2, num_G2, S_G2);
num_G3 = 1:length(fatigueValAfterORMSM_EMG_G3);
[p_G3, S_G3] = polyfit(num_G3, fatigueValAfterORMSM_EMG_G3, 2);
[fatigueValAfterORMNormSMFit_G3, delta_G3] = polyval(p_G3, num_G3, S_G3);
subplot(235), hold on
plot(fatigueValAfterORMNormSMFit_G1, 'b', 'LineWidth', 2)
plot(fatigueValAfterORMNormSMFit_G2, 'r', 'LineWidth', 2)
plot(fatigueValAfterORMNormSMFit_G3, 'm', 'LineWidth', 2)
grid on
title('2nd order fitting on EMG results')
xlabel('Workout #')
ylabel('Fatigue Level')

%% 3rd order Curve Fitting
num_G1 = 1:length(fatigueValAfterORMSM_G1);
[p_G1, S_G1] = polyfit(num_G1, fatigueValAfterORMSM_G1, 3);
[fatigueValAfterORMNormSMFit_G1, delta_G1] = polyval(p_G1, num_G1, S_G1);
num_G2 = 1:length(fatigueValAfterORMSM_G2);
[p_G2, S_G2] = polyfit(num_G2, fatigueValAfterORMSM_G2, 3);
[fatigueValAfterORMNormSMFit_G2, delta_G2] = polyval(p_G2, num_G2, S_G2);
num_G3 = 1:length(fatigueValAfterORMSM_G3);
[p_G3, S_G3] = polyfit(num_G3, fatigueValAfterORMSM_G3, 3);
[fatigueValAfterORMNormSMFit_G3, delta_G3] = polyval(p_G3, num_G3, S_G3);
subplot(233), hold on
plot(fatigueValAfterORMNormSMFit_G1, 'b', 'LineWidth', 2)
plot(fatigueValAfterORMNormSMFit_G2, 'r', 'LineWidth', 2)
plot(fatigueValAfterORMNormSMFit_G3, 'm', 'LineWidth', 2)
grid on
title('3rd order fitting on acoustic results')
xlabel('Workout #')
ylabel('Fatigue Level')


num_G1 = 1:length(fatigueValAfterORMSM_EMG_G1);
[p_G1, S_G1] = polyfit(num_G1, fatigueValAfterORMSM_EMG_G1, 3);
[fatigueValAfterORMNormSMFit_G1, delta_G1] = polyval(p_G1, num_G1, S_G1);
num_G2 = 1:length(fatigueValAfterORMSM_EMG_G2);
[p_G2, S_G2] = polyfit(num_G2, fatigueValAfterORMSM_EMG_G2, 3);
[fatigueValAfterORMNormSMFit_G2, delta_G2] = polyval(p_G2, num_G2, S_G2);
num_G3 = 1:length(fatigueValAfterORMSM_EMG_G3);
[p_G3, S_G3] = polyfit(num_G3, fatigueValAfterORMSM_EMG_G3, 3);
[fatigueValAfterORMNormSMFit_G3, delta_G3] = polyval(p_G3, num_G3, S_G3);
subplot(236), hold on
plot(fatigueValAfterORMNormSMFit_G1, 'b', 'LineWidth', 2)
plot(fatigueValAfterORMNormSMFit_G2, 'r', 'LineWidth', 2)
plot(fatigueValAfterORMNormSMFit_G3, 'm', 'LineWidth', 2)
grid on
title('3rd order fitting on EMG results')
xlabel('Workout #')
ylabel('Fatigue Level')

suptitle('Xingzhe on 07262019 with left hand 19 lbs.') % fill out the name*, data*, hand* and weight*
