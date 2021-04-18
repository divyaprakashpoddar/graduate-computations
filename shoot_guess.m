function x = shoot_guess(odefun,tspan,init,guess)
    [t,y] = ode45(odefun,tspan,[init guess]);
    x = y(end,1);
end
