clear; clc; close all;

data = csvread('Admission_Predict.csv',1,0);

x = data(:,2);
y = data(:,9);
n = numel(x);

scatter(x,y);
hold on

x2 = x.^2;
xy = x.*y;

% y = a + bx
a = (sum(y)*sum(x2)-sum(x)*sum(xy))/(n*sum(x2)-(sum(x))^2);
b = (n*sum(xy)-sum(x)*sum(y))/(n*sum(x2) - (sum(x))^2); 

xq = linspace(min(x),max(x),1000);
%plot(xq,a+b*xq,'r-.')

p = polyfit(x,y,1);
plot(xq,polyval(p,xq),'--');
