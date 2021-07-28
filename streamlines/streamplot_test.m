clear
close all
clc

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');


%% settings
Vinf = 1;
alpha = deg2rad(5);
c = 1;
N=100;
Nstream = 5;
streamRes = 1000;
xmin = -0.25;
xmax = 1.25;
ymin = -0.5;
ymax = 0.5;

%% solve flow
coords = genNACA4([2 4 12], c, N+1);
xc = coords(1:end-1,1) + diff(coords(:,1));
yc = coords(1:end-1,2) + diff(coords(:,2));


[q, g, ~, Cp] = panelHS(coords, Vinf, alpha);

%% plot coordinate generation

a = c/2;
Xe = xc - ones(size(xc))*c/2;
Ye = yc;

opts = optimset('Diagnostics','off', 'Display','off');
for i = 1:N
    [xsoln] = fsolve(@(X) ellip2cart([Xe(i),Ye(i)],X,a), [1,1], opts);
    mu(i) = xsoln(1);
    nu(i) = xsoln(2);
end

dispN = 30;
muDisp = [0, linspace(0.001,3,dispN-1)];
for i = 1:N
    for j = 1:dispN
        muD(i,j) = mu(i) + muDisp(j)*sign(mu(i));
        nuD(i,j) = nu(i);
    end
end
% muD = reshape(muD, 1, N*dispN);
% nuD = reshape(nuD, 1, N*dispN);
Xp = a.*cosh(muD).*cos(nuD) + ones(size(xc))*c/2;
Yp = a.*sinh(muD).*sin(nuD);

plotCoords = [reshape(Xp,N*dispN,1), reshape(Yp,N*dispN,1)];

%% calculate streamlines

[U, V, beta, r] = vField(coords, plotCoords, q, 1*g, Vinf, alpha);
Up = reshape(U, N, dispN);
Vp = reshape(V, N, dispN);
Up(:,1) = zeros(size(N,1)).*NaN;
Vp(:,1) = zeros(size(N,1)).*NaN;

[Xs, Ys] = meshgrid(linspace(xmin, xmin, Nstream), linspace(ymin, ymax, Nstream));
% Xs = rand(1,Nstream).*(xmax-xmin).^1 - ones(1,Nstream).*(xmax-xmin)/4;
% Ys = rand(1,Nstream).*(ymax-ymin).^1 - ones(1,Nstream).*(ymax-ymin)/2;


[xq, yq] = meshgrid(linspace(1.5*xmin,1.5*xmax,streamRes), linspace(1.5*ymin,1.5*ymax,streamRes));
uq = griddata(Xp*cos(-alpha)-Yp*sin(-alpha), Xp*sin(-alpha)+Yp*cos(-alpha), Up*cos(-alpha)-Vp*sin(-alpha), xq, yq);
vq = griddata(Xp*cos(-alpha)-Yp*sin(-alpha), Xp*sin(-alpha)+Yp*cos(-alpha), Up*sin(-alpha)+Vp*cos(-alpha), xq, yq);

%% plot streamlines and surface
figure('name', 'streamPlot', 'units', 'inches', 'papersize', [8 5], 'paperposition', [0 0 8 5])
ah = axes('Units','Normalize','Position',[0 0 1 1]);
hold on
AFpoly = polyshape([xc*cos(-alpha)-yc*sin(-alpha)], [xc*sin(-alpha)+yc*cos(-alpha)]);
af = plot(AFpoly);
af.FaceColor = [0.7 0.7 0.7];
af.FaceAlpha = 1;

axis equal
xlim([xmin,xmax])
ylim([ymin,ymax])
xlabel('x/c')
ylabel('z/c')
title('Free Flow')

% Cpplot = pcolor(xq, yq, ones(size(xq))-(uq.^2+vq.^2));
% Cpplot.FaceColor = 'interp';
% set(Cpplot, 'edgecolor', 'none')

%colorbar
caxis([0, 2])
[cpp, h] = contourf(xq, yq, sqrt(uq.^2+vq.^2), 10, 'color', 'none');
stream = streamline(stream2(xq, yq, uq, vq, [-0.25, -0.25], [-0.05, -0.10], [1, 1000]));
set(stream, 'color', 'black')
set(stream, 'linestyle', '-')
set(stream, 'linewidth', 2)

set(gca, 'Visible', 'off')

% view(2);

% shading interp
%print('Uniform3412.pdf', '-dpdf', '-painters')
print('velImage.png', '-dpng', '-r300')

%% functions
function F = ellip2cart(x,y,a)
F(1) = a*cosh(y(1))*cos(y(2)) - x(1);
F(2) = a*sinh(y(1))*sin(y(2)) - x(2);
end  
