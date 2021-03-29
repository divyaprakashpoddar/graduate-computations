function [c,iter] = secant(f,a,b,tol)
	iter = 0;
	c = (a*f(b) - b*f(a))/(f(b)-f(a));
	while abs(f(c)) > tol
		a = b;
		b = c;
		c = (a*f(b) - b*f(a))/(f(b)-f(a));
		iter = iter + 1;
	end
end
