clear
close all
clc

%% Parameters

alpha = 0;

Ns = [floor(logspace(1.5, 3, 10)/2)*2, 3000];


%% Constant Strength (midpoints)

for k = 1:length(Ns)
    
    V = 1;
    Qinf = [V*cos(alpha); V*sin(alpha)];

    N = Ns(k);
    coords = genNACA4([2 4 12], 1, N+1);

    xe = coords(:,1);
    ye = coords(:,2);
    xp = xe(2:end) - diff(xe)/2;
    yp = ye(2:end) - diff(ye)/2;
    thetap = atan2(diff(ye), diff(xe));
    np = [-sin(thetap), cos(thetap)];
    tp = [cos(thetap), sin(thetap)];
    ds = sqrt(diff(xp).^2 + diff(yp).^2);

    A = zeros(N+1);
    for i = 1:N
        for j = 1:N
            if j ~= i
                A(i, j) = constPhi(xp(i), yp(i), xe(j), ye(j), xe(j+1), ye(j+1), 1);
            else
                A(i, j) = 0.5;
            end
        end
        A(i, N+1) = constPhi(xp(i), yp(i), 1, 0, 1000, 0, 1);
    end

    A(N+1, :) = [1, zeros(1,N-2), -1, 1];

    RHS = -[[xp, yp]*Qinf; 0];

    mu = A\RHS;

    Qt = diff(mu(1:end-1))./ds;
    Cp = 1 - Qt.^2/(Qinf'*Qinf);
    
    Cl0(k) = -2*mu(end);
end

%% Continuous Linear Strength (midpoints)

for k = 1:length(Ns)
    V = 1;

    Qinf = [V*cos(alpha); V*sin(alpha)];

    N = Ns(k);
    coords = genNACA4([2 4 12], 1, N);

    xp = coords(:,1);
    yp = coords(:,2);

    s = [0; cumsum(sqrt(diff(xp).^2 + diff(yp).^2))];

    thetap = atan2(diff(-yp),-diff(xp));

    xc = xp(1:N-1) + 0.5*diff(xp);
    yc = yp(1:N-1) + 0.5*diff(yp);

    xc = xc;
    yc = yc;

    xc(end+1) = 1-(1e-9);
    yc(end+1) = 0;

    %

    Phia = zeros(N);
    Phib = zeros(N);

    % N panels, N+1 vertices (trailing edge doublet strength is not continuous)
    % collocation points i, panels j
    for i = 1:N
        for j = 1:N-1
            [Phia(i, j), Phib(i, j)] = linPhi(xc(i), yc(i), xp(j), yp(j), xp(j+1), yp(j+1), 1, 1);
        end
    end

    for i = 1:N-1
        Phia(i, i) = 0.25;
        Phib(i, i) = 0.25;
    end

    % assemble A
    A = Phia + circshift(Phib, 1, 2);

    % for i = 1:N
    %     A(i, i) = -0.25;
    % end

    % wake
    for i = 1:N
        A(i, N+1) = constPhi(xc(i), yc(i), 1, 0, 1000, 0, 1);
    end

    % kutta condition (TE potential jump applies to the wake)
    A(N+1, :) = [1, zeros(1, N-2), -1, 1];


    RHS = -[[xc, yc]*Qinf; 0];

    % kutta condition (TE doublet gradient constant)
    RHS(N) = 0;
    A(N, :) = [-1, 1, zeros(1, N-4), -1, 1, 0];

    mu = A\RHS;

    % take second order central differences to smooth oscillations
    Qt = gradient(mu(1:N))./gradient(s);
    Cp = 1 - Qt.^2/(Qinf'*Qinf);
    
    Cl1(k) = -2*mu(end);

end

%% Continuous Linear Strength (endpoints)

for k = 1:length(Ns)
    V = 1;

    Qinf = [V*cos(alpha); V*sin(alpha)];

    N = Ns(k);
    coords = genNACA4([2 4 12], 1, N);

    xp = coords(:,1);
    yp = coords(:,2);

    s = [0; cumsum(sqrt(diff(xp).^2 + diff(yp).^2))];

    thetap = atan2(diff(-yp),-diff(xp));

    thetac = [atan2(0.5*sin(thetap(1))+0.5*sin(thetap(end)), 0.5*cos(thetap(1))+0.5*cos(thetap(end))); ...
              atan2(0.5*sin(thetap(1:end-1))+0.5*sin(thetap(2:end)), 0.5*cos(thetap(1:end-1))+0.5*cos(thetap(2:end))); ...
              atan2(0.5*sin(thetap(1))+0.5*sin(thetap(end)), 0.5*cos(thetap(1))+0.5*cos(thetap(end)))];

    xc = xp - 1e-8*cos(thetac - pi/2);
    yc = yp - 1e-8*sin(thetac - pi/2);

    %

    Phia = zeros(N);
    Phib = zeros(N);

    % N panels, N+1 vertices (trailing edge doublet strength is not continuous)
    % collocation points i, panels j
    for i = 1:N
        for j = 1:N-1
            [Phia(i, j), Phib(i, j)] = linPhi(xc(i), yc(i), xp(j), yp(j), xp(j+1), yp(j+1), 1, 1);
        end
    end

    % assemble A
    A = Phia + circshift(Phib, 1, 2);


    % wake
    for i = 1:N
        A(i, N+1) = constPhi(xc(i), yc(i), 1, 0, 1000, 0, 1);
    end

    % kutta condition (TE potential jump applies to the wake)
    A(N+1, :) = [1, zeros(1, N-2), -1, 1];

    RHS = -[[xc, yc]*Qinf; 0];

    % kutta condition (TE doublet gradient constant)
    RHS(N) = 0;
    A(N, :) = [-1, 1, zeros(1, N-4), -1, 1, 0];

    mu = A\RHS;

    % take second order central differences to smooth oscillations
    Qt = diff(mu(1:N))./diff(s);
    Cp = 1 - Qt.^2/(Qinf'*Qinf);
    
    Cl1e(k) = -2*mu(end);

end

%% Discontinuous Linear Strength (gauss points)

for k = 1:length(Ns)
    V = 1;

    Qinf = [V*cos(alpha); V*sin(alpha)];

    N = Ns(k)/2;
    coords = genNACA4([2 4 12], 1, N+1);

    xp = coords(:,1);
    yp = coords(:,2);

    thetap = atan2(diff(-yp),-diff(xp));
    ds = sqrt(diff(xp).^2 + diff(yp).^2);

    xc1 = xp(1:N) + diff(xp)*(3-sqrt(3))/6; % - 1e-9*cos(thetap - pi/2);
    xc2 = xp(1:N) + diff(xp)*(3+sqrt(3))/6; % - 1e-9*cos(thetap - pi/2);
    yc1 = yp(1:N) + diff(yp)*(3-sqrt(3))/6; % - 1e-9*sin(thetap - pi/2);
    yc2 = yp(1:N) + diff(yp)*(3+sqrt(3))/6; % - 1e-9*sin(thetap - pi/2);

    xc = [xc1; xc2];
    yc = [yc1; yc2];

    s = [0; cumsum(sqrt(diff(xp).^2 + diff(yp).^2))];

    % figure
    % hold on
    % axis equal
    % grid on
    % grid minor
    % plot(xp, yp, 'marker', '.')
    % plot(xc1, yc1, 'linestyle', 'none', 'marker', '.')
    % plot(xc2, yc2, 'linestyle', 'none', 'marker', '.')

    % influences

    % first block: constant strengths
    % second block: linear strengths
    % points i, singularities j

    A = zeros(2*N+1);

    %      -----------------------------        |mu0(1)|
    %     | left points  | left points  |       | ...  |
    %     | constant     | linear       |       |mu0(N)|
    % A = |--------------+---------------  mu = |mu1(1)|
    %     | right points | right points |       | ...  |
    %     | constant     | linear       |       |mu1(N)|
    %      -----------------------------        |mu0(W)|

    % AICs
    for i = 1:N
        for j = 1:N
            if i == j
                A(i, j) = 0.5;
                A(i+N, j) = 0.5;
                A(i, j+N) = -sqrt(3)/6;
                A(i+N, j+N) = sqrt(3)/6;
            else
                % constant strengths
                A(i, j) = polyPhi(xc1(i), yc1(i), xp(j), yp(j), xp(j+1), yp(j+1), 1, 0, 0);
                A(i+N, j) = polyPhi(xc2(i), yc2(i), xp(j), yp(j), xp(j+1), yp(j+1), 1, 0, 0);
                % linear strengths
                A(i, j+N) = polyPhi(xc1(i), yc1(i), xp(j), yp(j), xp(j+1), yp(j+1), 0, 1, 0);
                A(i+N, j+N) = polyPhi(xc2(i), yc2(i), xp(j), yp(j), xp(j+1), yp(j+1), 0, 1, 0);
            end
        end
    end


    % wake
    for i = 1:N
        A(i, 2*N+1) = polyPhi(xc1(i), yc1(i), 1, 0, 1000, 0, 1, 0, 0);
        A(i+N, 2*N+1) = polyPhi(xc2(i), yc2(i), 1, 0, 1000, 0, 1, 0, 0);
    end

    % wake strength
    A(2*N+1, 1) = 1; A(2*N+1, N) = -1; A(2*N+1, N+1) = 0; A(2*N+1, 2*N) = 0;
    A(2*N+1, 2*N+1) = 1;

    % RHS
    RHS = -[[xc, yc]*Qinf; 0];

    % kutta velocity
%     RHS(2*N) = 0;
%     A(2*N, :) = zeros(1,2*N+1);
%     A(2*N, N+1) = ds(end); A(2*N, 2*N) = ds(1);


    % solve
    mu = A\RHS;

    % plot(xp(1:N) + diff(xp), mu(1:N))
    % 
    Qt = 2*mu(N+1:2*N)./ds;
    Cp = 1 - Qt.^2;

    Cl01(k) = -2*mu(end);
end

%% Convergence

figure
hold on
grid on
plot(Ns, abs([Cl0-Cl0(end); Cl1-Cl1(end); Cl1e-Cl1e(end); Cl01-Cl01(end)]))
plot(Ns, 1*Ns.^(-1), 'r:')
% plot(Ns, Ns.^(-1.5), 'k--')
plot(Ns, 30*Ns.^(-2), 'b:')
set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')

