clear
close all
clc

%% derivation

global N a b be c dt

N = 3;

[c, b] = lgwt(N, 0, 1); b = b';

C = zeros(N); V = zeros(N);
for i = 1:N
    for j = 1:N
        V(i, j) = c(i)^(j-1);
        C(i, j) = 1/j * c(i)^j;
    end
end
e = [1./(1:N-1), 0];

a = C/V;
be = e/V;

% main rule order 2N for gauss, 2N-1 for radau, 2N-2 for lobatto
% embedded error estimator order N-1

%% Solution

dt = 0.2;

t0 = 0;
t1 = 500;

nt = floor(t1/dt);

x0 = [1; 0; 0; 1.3];

t = [t0, dt*(1:nt)];
x = [x0, zeros(4, nt)];

ti = zeros(1, (N+1)*nt);
xi = zeros(4, (N+1)*nt);

opts = optimoptions('fsolve', 'SpecifyObjectiveGradient', false, 'Display','off', 'FunctionTolerance', 1e-15);

ks = zeros(4*N,1);
tic
for i = 1:nt
    %% Built-in Solver
    
%     ks = fsolve(@(k) kf(k, t(i), x(:,i)), [ks], opts);
    
    %% Newton Solver (faster, but with a timestep restriction)
    
    % euler step
    ke = kepler2d(t(i), x(:,i));
    
    % interpolate euler step to populate the initial guess of the stages
    for kj = 1:N
        ks(4*kj-(3:-1:0), 1) = kepler2d(t(i)+c(kj)*dt, x(:,i) + c(kj)*dt*ke);
    end
    
    % newton error vector
    Ner = kf(ks, t(i), x(:,i));
    
    % initialize
    iters = 0;
    J = zeros(4*N);
    j = zeros(4, 4, N);
    while norm(Ner) > 1e-14
        iters = iters+1;
        
        % build newton jacobian
        for Ji = 1:N
            krs = reshape(ks, [4, N])';
            [~, j(:,:,Ji)] = kepler2d(0, x(:,i) + a(N,:)*krs*dt);
            J(4*Ji-(3:-1:0),:) = dt*kron(a(Ji,:), j(:,:,Ji));
        end
        J = eye(4*N) - J;
        
        % next newton iteration
        k_next = ks - J\Ner;
        ks = k_next;
        Ner = kf(ks, t(i), x(:,i));
    end
    
    %% Stepper
    
    % take a step
    x(:,i+1) = x(:,i)' + dt*b*reshape(ks, [4 N])';
    
    % step error estimate
    Ser = dt*(b-be)*reshape(ks, [4 N])';
    
    % interpolate stages
    ti((N+1)*i-N) = t(i);
    xi(:,(N+1)*i-N) = x(:,i);

    ti((N+1)*i-(N-1:-1:0)) = t(i) + c*dt;
    xi(:,(N+1)*i-(N-1:-1:0)) =  x(:,i) + dt*reshape(ks, [4, N])*a';
end
ti(end+1) = t(end);
xi(:,end+1) = x(:,end);
GLNtime = toc;

Ei = 0.5*(xi(3,:).^2 + (xi(1,:).*xi(4,:)).^2) - 1./xi(1,:);

opts2 = odeset('reltol', 3e-14);
tic
[t45, x45] = ode45(@(t,x) kepler2d(t,x)', ti, x0, opts2);
RK45time = toc;
E45 = 0.5*(x45(:,3).^2 + (x45(:,1).*x45(:,4)).^2) - 1./x45(:,1);


figure
hold on
grid on
grid minor
plot(ti, Ei)
plot(t45, E45)
legend("GL"+2*N, "ODE45", 'location', 'best')
xlabel('Time')
ylabel('Total Energy')

abserr = norm(x45(end,:) - x(:,end)')

GLNtime/RK45time

%% functions

function [err, Jerr] = kf(k, t, x)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
global N a c dt

k = reshape(k, [4, N])';
err = zeros(4*N,1);
for i = 1:N
    err(4*i-(3:-1:0), 1) = k(i,:) - kepler2d(0, x' + dt*(a(i,:)*k));
end

if nargout == 2
    for Ji = 1:N
        [~, j(:,:,Ji)] = kepler2d(0, x' + a(N,:)*k*dt);
        J(4*Ji-(3:-1:0),:) = dt*kron(a(Ji,:), j(:,:,Ji));
    end
    Jerr = eye(4*N) - J;
end

end



function [xdot, J] = kepler2d(t, x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% x = [r, theta, rdot, thetadot];
% xdot = [rdot; thetadot; r*thetadot^2 - 1/r^2; -2/r * rdot * thetadot];

xdot = [x(3), x(4), x(1).*x(4).^2 - 1./(x(1).^2), -2./x(1) .* x(3) .* x(4)];

    if nargout == 2
    %     J =
    %     [                     0, 0,               1,            0]
    %     [                     0, 0,               0,            1]
    %     [    2/r^3 + thetadot^2, 0,               0, 2*r*thetadot]
    %     [ (2*rdot*thetadot)/r^2, 0, -(2*thetadot)/r,  -(2*rdot)/r]
        J = [0, 0, 1, 0; ...
             0, 0, 0, 1; ...
             2./(x(1).^3) + x(4).^2, 0, 0, 2.*x(1).*x(4); ...
            (2.*x(3).*x(4))./(x(1).^2), 0, -(2*x(4))./x(1), -(2*x(3))./x(1)];
    end
end

