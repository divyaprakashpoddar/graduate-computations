clear; clc; close all;

dx = 0.001;
x = 1:dx:2;
nx = numel(x);

% Boundary conditions
yI = 0;
yEnd = 0.5

odefun = @(x,y) [y(2); -exp(-x.*y(1)) - sin(y(2))];
g = fzero(@(guess) shoot_guess(odefun,x,yI,guess)-yEnd, 0.1);

% ODEs
dy1dx = @(x,y1,y2) y2;
dy2dx = @(x,y1,y2) -exp(-x.*y1) - sin(y2);

% Initialize
%guess = 0.1;
y1 = zeros(1,nx);
y2 = zeros(1,nx);
y1(1) = yI;
y2(1) = g;

for i = 1:nx-1
    y1(i+1) = y1(i) + dy1dx(x(i),y1(i),y2(i))*dx;
    y2(i+1) = y2(i) + dy2dx(x(i),y1(i),y2(i))*dx;
end

plot(x,y1)
hold on
plot(x(end),yEnd,'o')
