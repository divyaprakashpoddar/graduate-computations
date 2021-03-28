clear; clc; close all;

f = @(x) 4*log(x) - x;

xq = 0:0.1:20;
plot(xq,f(xq));
set(gca,'xAxisLocation','origin')

a = 1; b = 6;

c = (a*f(b) - b*f(a))/(f(b)-f(a));

tol = 1e-6;
counter = 0;
while abs(f(c)) > tol
	a = b;
	b = c;

	c = (a*f(b) - b*f(a))/(f(b)-f(a));
	counter = counter + 1;
end

