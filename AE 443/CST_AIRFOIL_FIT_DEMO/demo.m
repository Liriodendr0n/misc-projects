clc
clear all
close all

load RAE2822_FIT

Nfit = 20;

Winit = [-ones(1,Nfit), ones(1,Nfit)]; % initial weights

% Run the optimization code
[Wopt, fval]=fmincon(@(W) airfoilfit(W,yt,XL,XU,0),Winit,[],[],[],[],ones(1,2*Nfit)*-1,ones(1,2*Nfit),[]);

% Generate the CST airfoil with optimum weights
[ycoord] = CST_airfoil_fit(Wopt,XL,XU,0);

log10(fval)

% Plot and compare
plot(xcoord,ycoord,'b--')
hold on
plot(xcoord,yt,'r')
legend('CST','Target')
set(gcf,'color','w')
axis equal
