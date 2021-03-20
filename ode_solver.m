clear; clc; close all;

dydt = @(t) 2*t;

y0 = 1;

[t, y] = ode45(dydt,[0 5],y0);

plot(t,y,'-o','DisplayName','ode45')

yan = @(t) t.^2 + 1;
t = 0:0.1:5;

hold on
plot(t,yan(t),'--','linewidth',2,'DisplayName','Analytical')
xlabel('t')
ylabel('y')
legend('location','northwest')
