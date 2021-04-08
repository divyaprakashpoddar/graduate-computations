clear; clc; close all;

u = @(x,y,z,t) 2*x.^2 - x.*y*t.^2 + z.^2;
v = @(x,y,z,t) x.^2 - 4*x.*y + y.^2*t.^3;
w = @(x,y,z,t) 2*x.*y*t - y.*z + y.^2;

x = linspace(-1,1,8);
y = linspace(-1,1,8);
z = linspace(-1,1,8);

[x,y,z] = meshgrid(x,y,z);

for t = 0:0.1:5
    U = u(x,y,z,t);
    V = v(x,y,z,t);
    W = w(x,y,z,t);
    Umag = sqrt(U.^2+V.^2+W.^2);

    quiver3(x,y,z,U./Umag,V./Umag,W./Umag)
    drawnow
end
 
