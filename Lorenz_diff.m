clear; clc; close all;

rho = 28;
sigma = 10;
beta = 8/3;

odeFun = @(t,x) [sigma*(x(2)-x(1)); x(1)*(rho-x(3))-x(2); x(1)*x(2)-beta*x(3)];

trange = 0:0.01:100;
init = [1 1 1];

[t,x] = ode45(odeFun, trange, init);
plot3(x(:,1),x(:,2),x(:,3))

hold on
init = [1.000001 1 1];
[t,x] = ode45(odeFun, trange, init);
plot3(x(:,1),x(:,2),x(:,3))
