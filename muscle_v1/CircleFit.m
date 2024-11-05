%��дʱ�䣺2017.12.13
%��д��Ա��������
%��ϵ��ʽ��2220160354@bit.edu.cn
%ϵͳ�汾��Win10 ��ͥ���İ�
%����汾��R2015b
%�������ݣ�������С���˷���Բ����ϣ�����Ϊ�����꣬���ΪԲ������Ͱ뾶��
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

