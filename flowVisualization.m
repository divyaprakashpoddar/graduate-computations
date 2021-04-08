clear; clc; close all;

u = @(x,y,z,t) 2*x.^2*t - x.*y + z.^2;
v = @(x,y,z,t) x.^2 - 4*x.*y*t^2 + y.^2;
w = @(x,y,z,t) 2*x.*y - y.*z + y.^2*t^3;

x = linspace(-1,1,6);
y = linspace(-1,1,6);
z = linspace(-1,1,6);

[x,y,z] = meshgrid(x,y,z);
for t = 0:0.1:5
    U = u(x,y,z,t);
    V = v(x,y,z,t);
    W = w(x,y,z,t);

    Umag = sqrt(U.^2+V.^2+W.^2);

    quiver3(x,y,z,U./Umag,V./Umag,W./Umag)
    drawnow
end
