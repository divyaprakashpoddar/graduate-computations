clear; clc; close all;

f = @(x) 4*log(x) - x;

xq = 0:0.1:10;
plot(xq,f(xq));
set(gca,'xAxisLocation','origin')

a = 1; b = 6;

tol = 1e-6;

c = (a*f(b)-b*f(a))/(f(b)-f(a));

counter = 0;
while abs(f(c)) >  tol
	if f(a)*f(c) < 0
		b = c;
	else 
		a = c;
	end
	c = (a*f(b)-b*f(a))/(f(b)-f(a));
	counter = counter + 1;
end
