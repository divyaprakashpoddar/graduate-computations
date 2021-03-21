clear; clc; close all;

m = 10;
rho = 0.1;
g = 9.81;
A = 1;
cd = 0.5;

odeFun = @(t,y) [y(2); 1/m*(m*g - 1/2*rho*A*cd*y(2)^2)];
trange = [0 25];
H0 = 800; V0 = 0;
init = [H0 V0];

[t,y] = ode45(odeFun,trange,init)

figure(1)
plot(t,(H0 - (y(:,1)-H0)),'-o')
xlabel('time')
ylabel('position')
set(gca,'xaxislocation','origin')

figure(2)
plot(t,y(:,2),'-o')
xlabel('time')
ylabel('velocity')

