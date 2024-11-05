%编写时间：2017.12.13
%编写人员：李宇麒
%联系方式：2220160354@bit.edu.cn
%系统版本：Win10 家庭中文版
%软件版本：R2015b
%程序内容：基于最小二乘法的圆形拟合（输入为点坐标，输出为圆心坐标和半径）
%%
function [center_X, center_Y, R] = CircleFit(x, y)

%% Parameters
% the number of obtained points
Num = length(x);                                           
% defining parameters 
H =[x.',y.',ones(Num,1)];
Y = -1*[x.^2+y.^2].';

P = pinv(H'*H)*H'*Y;

center_X = -P(1)/2;
center_Y = -P(2)/2;
R = sqrt(center_X^2+center_Y^2-P(3));

