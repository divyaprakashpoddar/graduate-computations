clear; clc; close all;

u = @(x,y,z) 2*x.^2 - x.*y + z.^2;
v = @(x,y,z) x.^2 - 4*x.*y + y.^2;
w = @(x,y,z) 2*x.*y - y.*z + y.^2;

X = linspace(-1,1,6);
Y = linspace(-1,1,6);
Z = linspace(-1,1,6);

[x,y,z] = meshgrid(X,Y,Z);
U = u(x,y,z);
V = v(x,y,z);
W = w(x,y,z);

Umag = sqrt(U.^2+V.^2+W.^2);

quiver3(x,y,z,U./Umag,V./Umag,W./Umag)
xlabel('x'); ylabel('y'); zlabel('z')
hold on

xs = -0.5; ys = -0.5; zs = -0.5;

nSteps = 100;
dt = 0.1;
for t= 1:nSteps
    xs = xs + interpn(X,Y,Z,U,xs,ys,zs)*dt;
    ys = ys + interpn(X,Y,Z,V,xs,ys,zs)*dt;
    zs = zs + interpn(X,Y,Z,W,xs,ys,zs)*dt;
    if isnan(xs) || isnan(ys) || isnan(zs)
        break;
    end
    plot3(xs,ys,zs,'rx')
end
