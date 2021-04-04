clear; clc; close all;

EI = 100;
c = 1;
q = @(x) -c;
L = 1;

odefun = @(x,y) [y(2); y(3); y(4); q(x)/EI];

bcfun = @(ya,yb) [ya(1); yb(1); ya(3); yb(3)];

N = 10;
xmesh = linspace(0,L,N);
yGuess = ones(4,N);

solinit.x = xmesh;
solinit.y = yGuess;
  
sol = bvp4c(odefun,bcfun,solinit);

def = @(x) -c*x/24/EI.*(L^3-2*L*x.^2+x.^3);
slope = @(x) -c/24/EI.*(L^3-6*L*x.^2+4*x.^3);

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
