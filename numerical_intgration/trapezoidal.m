clear; clc; close all;

f = @(x) 2.^x;

% xq = 0:0.1:10;
% plot(xq,f(xq));

a = 0; b = 10;
N = 10000;
dx = (b-a)/N;

x = a:dx:b;
y = f(x);

tic
I = f(x(1)) + f(x(end));

for i = 2:N
    I = I + 2*f(x(i));
end

I = dx/2*I;
toc

tic
% Vectorization
IV = dx/2*(f(x(1)) + 2*sum(f(x(2:end-1))) + f(x(end)));
toc
