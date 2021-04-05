clear; clc; close all;

EI = 100;
F = 1;
q = @(x) 0;
L = 1;

odefun = @(x,y) [y(2); y(3); y(4); q(x)/EI];

bcfun = @(ya,yb) [ya(1); ya(2); yb(3); yb(4)-F/EI];

N = 10;
xmesh = linspace(0,L,N);
yGuess = ones(4,N);

solinit.x = xmesh;
solinit.y = yGuess;
  
sol = bvp4c(odefun,bcfun,solinit);

def = @(x) -F*x.^2/6/EI.*(3*L-x);
slope = @(x) -F*x/2/EI.*(2*L-x);

figure(1)
plot(sol.x,sol.y(1,:),'o','DisplayName','bvp4c')
hold on
plot(sol.x,def(sol.x),'DisplayName','Analytical')
legend show
title('Deflection')

figure(2)
plot(sol.x,sol.y(2,:),'o','DisplayName','bvp4c')
hold on
plot(sol.x,slope(sol.x),'DisplayName','Analytical')
legend show
title('Slope')
