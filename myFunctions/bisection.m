function [c, iter] = bisection(f,a,b,tol)
	iter = 0;
	c = mean([a,b]);
	while abs(f(c)) > tol
		[a,b] = rootInterval(f,a,b,c);
		c = mean([a,b]);
		iter = iter + 1;
	end
end
