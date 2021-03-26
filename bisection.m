clear; clc; close all;

f = @(x) 4*log(x) - x;

xq = 0:0.1:10;
plot(xq,f(xq))
set(gca,'xAxisLocation','origin')
hold on

a = 3; b = 9;

c = (a+b)/2;

tol = 1e-6;

while abs(f(c)) > tol
	if f(a)*f(c) < 0
		b = c;
	else
		a = c;
	end
	c = (a+b)/2;
	plot(c,0,'ro')
	pause(2)
end
plot(c,f(c),'rx')
