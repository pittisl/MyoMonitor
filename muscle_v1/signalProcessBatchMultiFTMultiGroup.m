% This is to plot multiple groups results 
%%%% This should process
% the data with MULTIPLE FEATURES
% THIS IS THE COPY OF signalProcessBatch.m
clear; clc;
fatigueValues_G1 = signalProcessBatchMultiFT;
% ----------------
fatigueValues_G2 = signalProcessBatchMultiFT;
fatigueValues_G3 = signalProcessBatchMultiFT;
% ----------------

%% Raw data plot
% figure, hold on 
% plot(fatigueValues_G1(1,:), 'b'), plot(fatigueValues_G1(2,:), '--b')
% plot(fatigueValues_G2(1,:), 'r'), plot(fatigueValues_G2(2,:), '--r')
% plot(fatigueValues_G3(1,:), 'm'), plot(fatigueValues_G3(2,:), '--m')
% grid on
% title('Raw fatigue measurements of two feature')
% xlabel('Workout #')
% ylabel('Fatigue Level')
% legend('Group 1 mean','Group 1 variance','Group 2 mean', 'Group 2 variance','Group 3 mean','Group 3 variance')

fatigueValues_G1_1 = fatigueValues_G1(1,:);
fatigueValues_G1_2 = fatigueValues_G1(2,:);
% ----------------
fatigueValues_G2_1 = fatigueValues_G2(1,:);
fatigueValues_G2_2 = fatigueValues_G2(2,:);
fatigueValues_G3_1 = fatigueValues_G3(1,:);
fatigueValues_G3_2 = fatigueValues_G3(2,:);
% ----------------

outlierIndex_G1_1 = isoutlier(fatigueValues_G1_1, 'median', 'ThresholdFactor', 6);
outlierIndex_G1_2 = isoutlier(fatigueValues_G1_2, 'median', 'ThresholdFactor', 6);
% ----------------
outlierIndex_G2_1 = isoutlier(fatigueValues_G2_1, 'median', 'ThresholdFactor', 6);
outlierIndex_G2_2 = isoutlier(fatigueValues_G2_2, 'median', 'ThresholdFactor', 6);
outlierIndex_G3_1 = isoutlier(fatigueValues_G3_1, 'median', 'ThresholdFactor', 6);
outlierIndex_G3_2 = isoutlier(fatigueValues_G3_2, 'median', 'ThresholdFactor', 6);
% ----------------
fatigueValAfterORM_G1_1 = fatigueValues_G1_1;
fatigueValAfterORM_G1_2 = fatigueValues_G1_2;
% ----------------
fatigueValAfterORM_G2_1 = fatigueValues_G2_1;
fatigueValAfterORM_G2_2 = fatigueValues_G2_2;
fatigueValAfterORM_G3_1 = fatigueValues_G3_1;
fatigueValAfterORM_G3_2 = fatigueValues_G3_2;
% ----------------
fatigueValAfterORM_G1_1(outlierIndex_G1_1) = [];
fatigueValAfterORM_G1_2(outlierIndex_G1_2) = [];
% ----------------
fatigueValAfterORM_G2_1(outlierIndex_G2_1) = [];
fatigueValAfterORM_G2_2(outlierIndex_G2_2) = [];
fatigueValAfterORM_G3_1(outlierIndex_G3_1) = [];
fatigueValAfterORM_G3_2(outlierIndex_G3_2) = [];
% ----------------
% figure, hold on 
% plot(fatigueValAfterORM_G1_1, 'b'), plot(fatigueValAfterORM_G1_2, '--b')
% plot(fatigueValAfterORM_G2_1, 'r'), plot(fatigueValAfterORM_G2_2, '--r')
% plot(fatigueValAfterORM_G3_1, 'm'), plot(fatigueValAfterORM_G3_2, '--m')
% grid on
% title('Raw fatigue measurements of two feature (after outlier removal)')
% xlabel('Workout #')
% ylabel('Fatigue Level')
% legend('Group 1 mean','Group 1 variance','Group 2 mean', 'Group 2 variance','Group 3 mean','Group 3 variance')

%% smooth
fatigueValAfterORMNormSM_G1_1 = smooth(fatigueValAfterORM_G1_1, 10)';
fatigueValAfterORMNormSM_G1_2 = smooth(fatigueValAfterORM_G1_2, 10)';
% ----------------
fatigueValAfterORMNormSM_G2_1 = smooth(fatigueValAfterORM_G2_1, 10)';
fatigueValAfterORMNormSM_G2_2 = smooth(fatigueValAfterORM_G2_2, 10)';
fatigueValAfterORMNormSM_G3_1 = smooth(fatigueValAfterORM_G3_1, 10)';
fatigueValAfterORMNormSM_G3_2 = smooth(fatigueValAfterORM_G3_2, 10)';
% ----------------
fatigueValAfterORMNormSM_G1_1(end) = [];fatigueValAfterORMNormSM_G1_1(1) = []; % remove the start and end point
fatigueValAfterORMNormSM_G1_2(end) = [];fatigueValAfterORMNormSM_G1_2(1) = [];
% ----------------
fatigueValAfterORMNormSM_G2_1(end) = [];fatigueValAfterORMNormSM_G2_1(1) = [];
fatigueValAfterORMNormSM_G2_2(end) = [];fatigueValAfterORMNormSM_G2_2(1) = [];
fatigueValAfterORMNormSM_G3_1(end) = [];fatigueValAfterORMNormSM_G3_1(1) = [];
fatigueValAfterORMNormSM_G3_2(end) = [];fatigueValAfterORMNormSM_G3_2(1) = [];
% ----------------
% figure, hold on 
% plot(fatigueValAfterORMNormSM_G1_1, 'b'), plot(fatigueValAfterORMNormSM_G1_2, '--b')
% plot(fatigueValAfterORMNormSM_G2_1, 'r'), plot(fatigueValAfterORMNormSM_G2_2, '--r')
% plot(fatigueValAfterORMNormSM_G3_1, 'm'), plot(fatigueValAfterORMNormSM_G3_2, '--m')
% grid on
% title('Raw fatigue measurements of two feature (after outlier removal and smoothing)')
% xlabel('Workout #')
% ylabel('Fatigue Level')
% legend('Group 1 mean','Group 1 variance','Group 2 mean', 'Group 2 variance','Group 3 mean','Group 3 variance')

% make the two features as the same length
fatigueValAfterORMNormSM_G1_1 = fatigueValAfterORMNormSM_G1_1(end-min([length(fatigueValAfterORMNormSM_G1_1), ... 
    length(fatigueValAfterORMNormSM_G1_2)]) + 1:end);
fatigueValAfterORMNormSM_G1_2 = fatigueValAfterORMNormSM_G1_2(end-min([length(fatigueValAfterORMNormSM_G1_1), ... 
    length(fatigueValAfterORMNormSM_G1_2)]) + 1:end);
% ----------------
fatigueValAfterORMNormSM_G2_1 = fatigueValAfterORMNormSM_G2_1(end-min([length(fatigueValAfterORMNormSM_G2_1), ... 
    length(fatigueValAfterORMNormSM_G2_2)]) + 1:end);
fatigueValAfterORMNormSM_G2_2 = fatigueValAfterORMNormSM_G2_2(end-min([length(fatigueValAfterORMNormSM_G2_1), ... 
    length(fatigueValAfterORMNormSM_G2_2)]) + 1:end);
fatigueValAfterORMNormSM_G3_1 = fatigueValAfterORMNormSM_G3_1(end-min([length(fatigueValAfterORMNormSM_G3_1), ... 
    length(fatigueValAfterORMNormSM_G3_2)]) + 1:end);
fatigueValAfterORMNormSM_G3_2 = fatigueValAfterORMNormSM_G3_2(end-min([length(fatigueValAfterORMNormSM_G3_1), ... 
    length(fatigueValAfterORMNormSM_G3_2)]) + 1:end);
% ----------------

%% linear regression
y = 1:length(fatigueValAfterORMNormSM_G1_1);
x1 = fatigueValAfterORMNormSM_G1_1';
x2 = fatigueValAfterORMNormSM_G1_2';

X = [ones(size(x1)) x1 x2 x1.^2 x2.^2 x1.*x2];
b1 = regress(y',X);

figure
subplot(131)
scatter3(x1,x2,y,'filled')
hold on
x1fit = min(x1):0.01:max(x1);
x2fit = min(x2):0.01:max(x2);
[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
YFIT = b1(1) + b1(2)*X1FIT + b1(3)*X2FIT + b1(4)*X1FIT.^2 + b1(5)*X2FIT.^2 + b1(6)*X1FIT.*X2FIT;
mesh(X1FIT,X2FIT,YFIT)
xlabel('Mean')
ylabel('Skewness')
zlabel('Workout Index #'),title('Group 1')
view(50,10)
% ----------------
%
y = 1:length(fatigueValAfterORMNormSM_G2_1);
x1 = fatigueValAfterORMNormSM_G2_1';
x2 = fatigueValAfterORMNormSM_G2_2';

X = [ones(size(x1)) x1 x2 x1.^2 x2.^2 x1.*x2];
b2 = regress(y',X);


% figure
subplot(132)
scatter3(x1,x2,y,'filled')
hold on
x1fit = min(x1):0.01:max(x1);
x2fit = min(x2):0.01:max(x2);
[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
YFIT = b2(1) + b2(2)*X1FIT + b2(3)*X2FIT + b2(4)*X1FIT.^2 + b2(5)*X2FIT.^2 + b2(6)*X1FIT.*X2FIT;
mesh(X1FIT,X2FIT,YFIT)
xlabel('Mean')
ylabel('Skewness')
zlabel('Workout Index #'), title('Group 2')
view(50,10)

%
y = 1:length(fatigueValAfterORMNormSM_G3_1);
x1 = fatigueValAfterORMNormSM_G3_1';
x2 = fatigueValAfterORMNormSM_G3_2';

X = [ones(size(x1)) x1 x2 x1.^2 x2.^2 x1.*x2];
b3 = regress(y',X);

% figure
subplot(133)
scatter3(x1,x2,y,'filled')
hold on
x1fit = min(x1):0.01:max(x1);
x2fit = min(x2):0.01:max(x2);
[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
YFIT = b3(1) + b3(2)*X1FIT + b3(3)*X2FIT + b3(4)*X1FIT.^2 + b3(5)*X2FIT.^2 + b3(6)*X1FIT.*X2FIT;
mesh(X1FIT,X2FIT,YFIT)
xlabel('Mean')
ylabel('Skewness')
zlabel('Workout Index #'),title('Group 3')
view(50,10)