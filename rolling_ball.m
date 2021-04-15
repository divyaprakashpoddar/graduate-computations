clear; clc; close all;

syms x y

f = 3*(1-x)^2*exp(-x^2 - (y+1)^2) ...
              - 10*(x/5 - x^3 - y^5)*exp(-x^2-y^2) ...
              - 1/3*exp(-(x+1)^2 - y^2);

Gx = gradient(f,x);
Gy = gradient(f,y);

f = matlabFunction(f);
Gx = matlabFunction(Gx);
Gy = matlabFunction(Gy);

N = 50;
x = linspace(-2.5,2.5,N);
y = linspace(-2.5,2.5,N);

[x,y] = meshgrid(x,y);

surf(x,y,f(x,y))
hold on

xpos = -0.5; ypos = 1.5;
plot3(xpos,ypos,f(xpos,ypos),'ro','markerfacecolor','r')

nSteps = 100;
dx = 0.01; dy = 0.01;
for i = 1:nSteps
    xTemp = xpos - Gx(xpos,ypos)*dx; 
    yTemp = ypos - Gy(xpos,ypos)*dy;

    xpos = xTemp; ypos = yTemp;
    plot3(xpos,ypos,f(xpos,ypos),'ro','markerfacecolor','r')
    pause(0.5)
end 
