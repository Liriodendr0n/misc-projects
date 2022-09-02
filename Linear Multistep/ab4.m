% Usage: [y t] = ab4(f,a,b,ya,n) or y = ab4(f,a,b,ya,n)
% Adams-Bashforth 4-step method for initial value problems
% It uses Runge-Kutta method of order 4 as a starter
%
% Input:
% f - Matlab inline function f(t,y)
% a,b - interval
% ya - initial condition
% n - number of subintervals (panels)
%
% Output:
% y - computed solution
% t - time steps
%
% Examples:
% [y t]=ab4(@myfunc,0,1,1,10);          here 'myfunc' is a user-defined function in M-file
% y=ab4(inline('sin(y*t)','t','y'),0,1,1,10);
% f=inline('sin(y(1))-cos(y(2))','t','y');
% y=ab4(f,0,1,1,10);

function [t y] = ab4(f,a,b,ya,n)
h = (b - a) / n;
h24 = h / 24;

y(1,:) = ya;
t(1) = a;

m = min(3,n);

for i = 1 : m % start-up phase, using Runge-Kutta of order 4
    t(i+1) = t(i) + h;
    s(i,:) = f(t(i), y(i,:));
    s2 = f(t(i) + h / 2, y(i,:) + s(i,:) * h /2);
    s3 = f(t(i) + h / 2, y(i,:) + s2 * h /2);
    s4 = f(t(i+1), y(i,:) + s3 * h);
    y(i+1,:) = y(i,:) + (s(i,:) + s2+s2 + s3+s3 + s4) * h / 6;
end;

for i = m + 1 : n % main phase
    s(i,:) = f(t(i), y(i,:));
    y(i+1,:) = y(i,:) + (55 * s(i,:) - 59 * s(i-1,:) + 37 * s(i-2,:) - 9 * s(i-3,:)) * h24;
    t(i+1) = t(i) + h;
end;

y = y';