%starting point
x = [ 10.5440; 4.1124; 35.8233];

dt = 0.01;
N=10000;
x_series = [x];
for i=1:N
  x = gauss_step( x, @lorenz_dynamics, dt, 1e-7, 1, 100);
  x_series = [x_series x];
end

plot3( x_series(1,:), x_series(2,:), x_series(3,:) );
set(gca,'xtick',[],'ytick',[],'ztick',[]);
title('Lorenz Attractor');
return;

function [td, j] = lorenz_dynamics(state)
  %return a time derivative and a Jacobian of that time derivative
  x = state(1);
  y = state(2);
  z = state(3);
  
  sigma = 10;
  beta  = 8/3;
  rho   = 28;
  
  td = [sigma*(y-x); x*(rho-z)-y; x*y-beta*z];

  j = [-sigma, sigma, 0;
        rho-z, -1, -x;
        y, x, -beta];
end

function x_next = gauss_step( x, dynamics, dt, threshold, damping, max_iterations )
  [d,~] = size(x);
  sq3 = sqrt(3);
  if damping > 1 || damping <= 0
    error('damping should be between 0 and 1.') 
  end
  
  %Use explicit Euler steps as initial guesses
  [k,~] = dynamics(x);
  x1_guess = x + (1/2-sq3/6)*dt*k;
  x2_guess = x + (1/2+sq3/6)*dt*k;
  [k1,~] = dynamics(x1_guess);
  [k2,~] = dynamics(x2_guess);
  
  a11 = 1/4;
  a12 = 1/4 - sq3/6;
  a21 = 1/4 + sq3/6;
  a22 = 1/4;
  
  error = @(k1,k2) [ k1 - dynamics(x+(a11*k1+a12*k2)*dt); k2 - dynamics(x+(a21*k1+a22*k2)*dt) ];
  er = error(k1,k2);
  iteration=1;
  while( norm(er) > threshold && iteration < max_iterations )
    %fprintf('Newton iteration %d: error is %f.\n', iteration, norm(er) );
    iteration = iteration + 1;
    
    [~, j1] = dynamics(x+(a11*k1+a12*k2)*dt);
    [~, j2] = dynamics(x+(a21*k1+a22*k2)*dt);
    
    j = [ eye(d) - dt*a11*j1, -dt*a12*j1;
         -dt*a21*j2, eye(d) - dt*a22*j2 ];
     
    k_next = [k1;k2] - damping*linsolve(j,er);
    k1 = k_next(1:d);
    k2 = k_next(d+(1:d));
    
    er = error(k1,k2);
  end
  if norm(er) > threshold 
    error('Newton did not converge by %d iterations.', max_iterations); 
  end
  x_next = x + dt/2*(k1+k2);
end