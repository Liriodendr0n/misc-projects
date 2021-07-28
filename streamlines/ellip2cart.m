function F = ellip2cart(x,y,a)

F(1) = a*cosh(y(1))*cos(y(2)) - x(1);
F(2) = a*sinh(y(1))*sin(y(2)) - x(2);

end

