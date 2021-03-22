clear; clc; close all;

m = 1;
k = 1;
c = input('Damping, C: ');

F = @(t) (t>=0 && t<=100)*5;

odeFun = @(t,y) [y(2); 1/m*(F(t)-c*y(2)-k*y(1))];
trange = 0:0.1:50;
y0 = 0; v0 = 0;
init = [y0 v0];

[t, y] = ode45(odeFun, trange, init);

figure(1)
plot(t,y(:,1),'-o')
xlabel('Time')
ylabel('Displacement')

figure(2)
plot(t,y(:,2),'-o')
xlabel('Time')
ylabel('Velocity')
