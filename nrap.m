clear; clc; close all;

syms x

f = 4*log(x) - x;
df = diff(f,x);

f = matlabFunction(f);
df = matlabFunction(df);

x = 1;

xq = 0:0.1:10;
plot(xq,f(xq))

%for i = 1:10
%	x = x - f(x)/df(x);
%end

counter = 0;
tol = 1e-6;
while abs(f(x)) > tol
	x = x - f(x)/df(x);
	counter = counter + 1;
end

hold on
plot(x,f(x),'rx')
set(gca,'xaxislocation','origin')
