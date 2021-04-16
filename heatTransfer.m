clear; clc; close all;

Nx = 10; Ny = 10;

T = zeros(Nx,Ny);

% Boundary conditions
T(1,:) = 100;
T(:,1) = 0;
T(:,end) = 0;
T(end,:) = 0;

TOld = T;

err = 1;
tol = 1e-10;

while err > tol
    for i = 2:Nx-1
        for j = 2:Ny-1
            TOld(i,j) = 0.25*(T(i+1,j) + T(i-1,j) + T(i,j+1) + T(i,j-1));
            err = norm(TOld - T);
            T = TOld;
        end
    end
end

contourf(T,'linecolor','none')
colorbar
hold on
[gx, gy] = gradient(T);
quiver(-gx,-gy)
