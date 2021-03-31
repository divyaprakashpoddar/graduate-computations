function I = my_trapz(f,a,b,N)
  dx = (b-a)/N;
  x = a:dx:b;
  I = dx/2*(f(x(1)) + 2*sum(f(x(2:end-1))) + f(x(end)));
end
