function [a,b] = rootInterval(f,a,b,c)
	if f(a)*f(c) < 0
		b = c;
	else
		a = c;
	end
end
