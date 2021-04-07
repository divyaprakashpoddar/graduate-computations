clear; clc; close all;

data = csvread('~/Downloads/position_salaries.csv',1,1);

x = data(:,1);
y = data(:,2);

k = 4;

V = x.^(0:k);

a = pinv(V)*y;

xq = linspace(x(1),x(end),100);

plot(x,y,'r*')
hold on
plot(xq,polyval(flip(a),xq))
