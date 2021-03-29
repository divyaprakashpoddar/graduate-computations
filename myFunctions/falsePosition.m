function [c,iter] = falsePosition(f,a,b,tol)
	iter = 0;
	c = (a*f(b) - b*f(a))/(f(b)-f(a));
	while abs(f(c)) > tol
		[a,b] = rootInterval(f,a,b,c);
		c = (a*f(b) - b*f(a))/(f(b)-f(a));
		iter = iter + 1;
	end
end
