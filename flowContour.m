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

%quiver3(x,y,z,U./Umag,V./Umag,W./Umag)
%xlabel('x'); ylabel('y'); zlabel('z')

% Contour plots
figure(1)
nz = 2;
contour(x(:,:,nz),y(:,:,nz),U(:,:,nz),'showtext','on',10)
figure(2)
contourf(x(:,:,nz),y(:,:,nz),Umag(:,:,nz))
colormap(jet)

% Slice 
figure(3)
slice(X,Y,Z,Umag,-0.5,-0.5,-0.5)
colormap(lines)
