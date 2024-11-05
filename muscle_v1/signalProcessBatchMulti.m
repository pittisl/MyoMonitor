% This is to plot multiple groups results
clear; clc;
fatigueValues_G1 = signalProcessBatch;
% ----------------
fatigueValues_G2 = signalProcessBatch;
fatigueValues_G3 = signalProcessBatch;
% ----------------

%% Raw data plot
% figure, hold on 
% plot(fatigueValues_G1, 'b')
% plot(fatigueValues_G2, 'r')
% plot(fatigueValues_G3, 'm')
% grid on
% title('Raw fatigue measurements')
% xlabel('Workout #')
% ylabel('Fatigue Level')
% legend('Group 1','Group 2','Group 3')

%% outlier
outlierIndex_G1 = isoutlier(fatigueValues_G1, 'median', 'ThresholdFactor', 6);
outlierIndex_G2 = isoutlier(fatigueValues_G2, 'median', 'ThresholdFactor', 6);
outlierIndex_G3 = isoutlier(fatigueValues_G3, 'median', 'ThresholdFactor', 6);
fatigueValAfterORM_G1 = fatigueValues_G1;
fatigueValAfterORM_G2 = fatigueValues_G2;
fatigueValAfterORM_G3 = fatigueValues_G3;
fatigueValAfterORM_G1(outlierIndex_G1) = [];
fatigueValAfterORM_G2(outlierIndex_G2) = [];
fatigueValAfterORM_G3(outlierIndex_G3) = [];
% figure, hold on 
% plot(fatigueValAfterORM_G1, '--ob')
% plot(fatigueValAfterORM_G2, '--or')
% plot(fatigueValAfterORM_G3, '--om')
% grid on
% title('Raw fatigue measurements without outliers')
% xlabel('Workout #')
% ylabel('Fatigue Level')

%% normalization
% fatigueValAfterORMNorm_G1 = mapminmax(fatigueValAfterORM_G1);
% fatigueValAfterORMNorm_G2 = mapminmax(fatigueValAfterORM_G2);
% fatigueValAfterORMNorm_G3 = mapminmax(fatigueValAfterORM_G3);
% figure, hold on 
% plot(fatigueValAfterORMNorm_G1, '--ob')
% plot(fatigueValAfterORMNorm_G2, '--or')
% plot(fatigueValAfterORMNorm_G3, '--om')
% grid on
% axis([-inf inf -2 2]);
% title('Normalized raw fatigue measurements without outliers')
% xlabel('Workout #')
% ylabel('Fatigue Level')

%% smooth
fatigueValAfterORMNormSM_G1 = smooth(fatigueValAfterORM_G1, 10)';
% ----------------
fatigueValAfterORMNormSM_G2 = smooth(fatigueValAfterORM_G2, 10)';
fatigueValAfterORMNormSM_G3 = smooth(fatigueValAfterORM_G3, 10)';
% ----------------
fatigueValAfterORMNormSM_G1(end) = [];fatigueValAfterORMNormSM_G1(1) = []; % remove the start and end point
% ----------------
fatigueValAfterORMNormSM_G2(end) = [];fatigueValAfterORMNormSM_G2(1) = [];
fatigueValAfterORMNormSM_G3(end) = [];fatigueValAfterORMNormSM_G3(1) = [];
% ----------------

% figure, hold on 
% plot(fatigueValAfterORMNormSM_G1, 'xb')
% plot(fatigueValAfterORMNormSM_G2, 'xr')
% plot(fatigueValAfterORMNormSM_G3, 'xm')
% grid on
% title('Normalized and Smoothed raw fatigue measurements without outliers')
% xlabel('Workout #')
% ylabel('Fatigue Level')

%% segmented curve fitting
% ----------------
curveSearch(fatigueValAfterORMNormSM_G1, fatigueValAfterORMNormSM_G2, fatigueValAfterORMNormSM_G3);
% ----------------

%% Curve Fitting
% num_G1 = 1:length(fatigueValAfterORMNormSM_G1);
% [p_G1, S_G1] = polyfit(num_G1, fatigueValAfterORMNormSM_G1, 3);
% [fatigueValAfterORMNormSMFit_G1, delta_G1] = polyval(p_G1, num_G1, S_G1);
% num_G2 = 1:length(fatigueValAfterORMNormSM_G2);
% [p_G2, S_G2] = polyfit(num_G2, fatigueValAfterORMNormSM_G2, 3);
% [fatigueValAfterORMNormSMFit_G2, delta_G2] = polyval(p_G2, num_G2, S_G2);
% num_G3 = 1:length(fatigueValAfterORMNormSM_G3);
% [p_G3, S_G3] = polyfit(num_G3, fatigueValAfterORMNormSM_G3, 3);
% [fatigueValAfterORMNormSMFit_G3, delta_G3] = polyval(p_G3, num_G3, S_G3);
% figure, hold on
% plot(fatigueValAfterORMNormSMFit_G1, 'b', 'LineWidth', 2)
% plot(fatigueValAfterORMNormSMFit_G2, 'r', 'LineWidth', 2)
% plot(fatigueValAfterORMNormSMFit_G3, 'm', 'LineWidth', 2)
% grid on
% title('Normalized and Smoothed raw fatigue measurements without outliers with fitting results')
% xlabel('Workout #')
% ylabel('Fatigue Level')

%% Nonlinear regression
% X_G1 = num_G1;
% Y = fatigueValAfterORMNormSM_G1;
% modelfun = @(b, x)(b(1).*x.^3 + b(2).*x.^2 + b(3).*x + b(4)); %+ b(4)*exp(b(5).*x));
% beta0 = [0.1 0.2 0.3 0.4 ];
% 
% opts = statset('nlinfit');
% opts.RobustWgtFun = 'bisquare';
% 
% beta_G1 = nlinfit(X_G1,Y,modelfun,beta0, opts);
% X_G2 = num_G2;
% Y = fatigueValAfterORMNormSM_G2;
% beta_G2 = nlinfit(X_G2,Y,modelfun,beta0, opts);
% X_G3 = num_G3;
% Y = fatigueValAfterORMNormSM_G3;
% beta_G3 = nlinfit(X_G3,Y,modelfun,beta0, opts);
% 
% 
% YY_G1 = modelfun(beta_G1, X_G1);
% YY_G2 = modelfun(beta_G2, X_G2);
% YY_G3 = modelfun(beta_G3, X_G3);
% 
% figure, hold on
% plot(YY_G1, 'b', 'LineWidth', 2)
% plot(YY_G2, 'r', 'LineWidth', 2)
% plot(YY_G3, 'm', 'LineWidth', 2)
% grid on

%
% figure, hold on
% plot(repreValueFatigue_wo_out_rev,'x' , 'LineWidth', 2.0)
% repreValueFatigue_wo_out_sm = smooth(repreValueFatigue_wo_out_rev)';
% plot(repreValueFatigue_wo_out_sm,'LineWidth', 2.0)
% num = 1:length(repreValueFatigue_wo_out_sm);
% [p, S] = polyfit(num, repreValueFatigue_wo_out_sm, 3);
% [repreValueFatigue_fit, delta] = polyval(p, num, S);
% plot(repreValueFatigue_fit, 'LineWidth', 2.0)
% legend('Raw','5-point average','curve fitting'),grid on
% % axis([0 inf 0 10]);
% divider_indices = strfind(filepath,'\');
% identifier = filepath(divider_indices(end)+1:end);
% title(['Fatigue Indicator @', identifier, ' (error = ', num2str(sum(abs(delta))) ,')']), 
% xlabel('Workout #')
% ylabel('Fatigue Level')
% % compute the slope
% slope = (repreValueFatigue_fit(end) - repreValueFatigue_fit(1))/(length(repreValueFatigue_fit) - 1);
% annotation('textbox', [.2 .2 .7 .7], 'String', ['Slope = ', num2str(slope)], 'FitBoxToText', 'on')