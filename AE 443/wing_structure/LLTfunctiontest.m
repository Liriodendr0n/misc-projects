clear
close all
clc

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% test parameters

rho = 0.0023;
CLopt = 0.6;
W = 14000;

wing.AR = 6.5;
wing.S = 200;
wing.lambda = 0.5;
wing.kappa = 0.01;

flight.v_inf = sqrt(2*W/(rho*wing.S*CLopt));
flight.rho_inf = 0.00238;
flight.alpha = deg2rad(4.7);
flight.omega = deg2rad(0);

section.a0Lflap = deg2rad(-3);
section.a0LailR = deg2rad(-3);
section.a0LailL = deg2rad(-3);

b = sqrt(wing.AR*wing.S);
c0 = 2*b/(wing.AR*(wing.kappa + wing.lambda - wing.kappa*wing.lambda + 1));

[Lprime, cs, ys, alpha0L] = LLT(wing, section, flight);

%% shear force and bending moment

rollM = trapz(ys, Lprime.*ys)

ShearF = -cumtrapz(sign(ys).*ys, Lprime);
BendM = cumtrapz(sign(ys).*ys, ShearF).*sign(ys);
%figure
%hold on
%plot(ys, Lprime)
%plot(ys, ShearF)
%plot(ys, BendM)

figure
plot(ys, Lprime./(0.5*flight.rho_inf*flight.v_inf^2.*cs))
ylim([0 1])
% plot(ys, 2*Lprime./(rho.*cs.*flight.v_inf.^2))
% plot(ys, Lprime)

L = trapz(ys, Lprime)
CL = 2*L/(flight.rho_inf*wing.S*flight.v_inf^2)
% 
% figure('units', 'inches', 'papersize', [7 3])
% 
% hold on
% grid on
% grid minor
% plot(ys, Lprime)
% xlim([-b/2 b/2])
% ylim([0 4000])
% ylabel('Lprime [lbf/ft]')
% xlabel('Span Coordinate [ft]')
% 
% print('Lprime_maxload.pdf', '-dpdf', '-painters', '-fillpage')
% 
% figure('units', 'inches', 'papersize', [7 3])
% hold on
% grid on
% grid minor
% plot(ys, ShearF)
% xlim([-b/2 b/2])
% ylim([0 55000])
% ylabel('Shear Force [lbf]')
% xlabel('Span Coordinate [ft]')
% ax = gca;
% ax.YAxis.Exponent = 0;
% 
% print('ShearF_maxload.pdf', '-dpdf', '-painters', '-fillpage')
% 
% figure('units', 'inches', 'papersize', [7 3])
% hold on
% grid on
% grid minor
% plot(ys, BendM)
% xlim([-b/2 b/2])
% ylim([-400000 400000])
% ylabel('Bending Moment [lbf ft]')
% xlabel('Span Coordinate [ft]')
% ax = gca;
% ax.YAxis.Exponent = 0;
% 
% print('BendM_maxload.pdf', '-dpdf', '-painters', '-fillpage')

% figure('units', 'inches', 'papersize', [4 4])
% hold on
% grid on
% grid minor
% 
% yyaxis left
% plot(ys, ShearF)
% xlim([-b/2 b/2])
% ylim([-60000 60000])
% ylabel('Shear Force [lbf]')
% ax = gca;
% ax.YAxis(1).Exponent = 0;
% 
% yyaxis right
% plot(ys, BendM)
% xlim([-b/2 b/2])
% ylim([-400000 400000])
% ylabel('Bending Moment [lbf ft]')
% 
% ax = gca;
% ax.YAxis(2).Exponent = 0;
% xlabel('Span Position [ft]')
% 
% legend('$V_{wing}$', '$M_{wing}$', 'location', 'southwest')
% 
% %print('shearbendPDP.png', '-dpng', '-r600')
% 
% %% visualization
% figure('units', 'inches', 'papersize', [4 2])
% hold on
% axis equal
% grid on
% view([1,-0.5,0.33])
% 
% % lift per span
% quiver3(zeros(size(ys)), ys, zeros(size(ys)), zeros(size(Lprime)), zeros(size(Lprime)), Lprime)
% %lifting line
% plot3([0, 0], [-b/2, b/2], [0, 0], 'color', 'black')
% % wing outline + flaps
% % LE
% plot3([-0.25*c0*wing.lambda, -0.25*c0, -0.25*c0, -0.25*c0*wing.lambda], [-b/2, -wing.kappa*b/2, wing.kappa*b/2, b/2], [0,0,0,0], 'color', 'black')
% % 0.75c
% plot3([0.50*c0*wing.lambda, 0.50*c0, 0.50*c0, 0.50*c0*wing.lambda], [-b/2, -wing.kappa*b/2, wing.kappa*b/2, b/2], [0,0,0,0], 'color', 'black')
% % left ail
% plot3([0.75*c0*wing.lambda, 0.75*c0], [-b/2, -wing.kappa*b/2], [wing.lambda*0.25*c0*section.a0LailL, 0.25*c0*section.a0LailL], 'color', 'black')
% % right ail
% plot3([0.75*c0, 0.75*c0*wing.lambda], [wing.kappa*b/2, b/2], [0.25*c0*section.a0LailR, wing.lambda*0.25*c0*section.a0LailR], 'color', 'black')
% % flap
% plot3([0.75*c0, 0.75*c0], [-wing.kappa*b/2, wing.kappa*b/2], [0.25*c0*section.a0Lflap, 0.25*c0*section.a0Lflap], 'color', 'black')
% % tips
% plot3(wing.lambda*[-0.25*c0, 0, 0.5*c0, 0.75*c0], -b/2*[1,1,1,1], wing.lambda*[0,0,0,0.25*c0*section.a0LailL], 'color', 'black')
% plot3(wing.lambda*[-0.25*c0, 0, 0.5*c0, 0.75*c0], b/2*[1,1,1,1], wing.lambda*[0,0,0,0.25*c0*section.a0LailR], 'color', 'black')
% % outer breaks
% plot3([-0.25*c0, 0, 0.5*c0, 0.75*c0], -wing.kappa*b/2*[1,1,1,1], [0,0,0,0.25*c0*section.a0LailL], 'color', 'black')
% plot3([-0.25*c0, 0, 0.5*c0, 0.75*c0], wing.kappa*b/2*[1,1,1,1], [0,0,0,0.25*c0*section.a0LailR], 'color', 'black')
% % inner breaks
% plot3([-0.25*c0, 0, 0.5*c0, 0.75*c0], -wing.kappa*b/2*[1,1,1,1], [0,0,0,0.25*c0*section.a0Lflap], 'color', 'black')
% plot3([-0.25*c0, 0, 0.5*c0, 0.75*c0], wing.kappa*b/2*[1,1,1,1], [0,0,0,0.25*c0*section.a0Lflap], 'color', 'black')
% 
% ylim([-b/2, b/2])
% xlim([-c0/2, 3*c0/2])
% zlim([-3, 7])
% 
% xlabel('x [ft]')
% ylabel('y [ft]')
% zlabel('z [ft]')
% 
% %print('LLTview_wingroll_neg120roll.pdf', '-dpdf', '-painters', '-fillpage')
% %print('wingloadPDP.png', '-dpng', '-r600')
%%


