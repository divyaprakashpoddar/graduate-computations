function I = integralSum(f,a,b,n,typeI)
	dx = (b-a)/n;
	x = a:dx:b;
	if strcmp(typeI,'left')
		% Left Riemann Sum
		I = sum(dx*f(x(1:n)));
	elseif strcmp(typeI,'right')
		% Right Riemann Sum
		I = sum(dx*f(x(2:n+1)));
	elseif strcmp(typeI,'mid')
		% Mid-Point Method
		c = (a+dx/2):dx:(b-dx/2);
		I = sum(dx*f(c));
	else
		printf('Invalid integration type!\n')
	end
end
