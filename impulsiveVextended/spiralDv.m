clear
close all
clc

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

alphas = logspace(0, -1, 300);

parfor i = 1:length(alphas)
    [Dvesc_tan(i), tesc_tan(i), gammaesc_tan(i), thetaesc_tan(i)] = escDv_tan(alphas(i));
    [Dvesc_azi(i), tesc_azi(i), gammaesc_azi(i), thetaesc_azi(i)] = escDv_azi(alphas(i));
    % i/length(alphas)
end

%% delta v
figure('name', 'v', 'units', 'inches', 'papersize', [6 6], 'paperposition', [0 0 6 6])
grid on
grid minor
hold on
plot(alphas, Dvesc_tan)
plot(alphas, Dvesc_azi)
ylim([0.4, 1])
xlim([alphas(end), alphas(1)])
set(gca, 'Xscale', 'log')
set(gca, 'Xdir', 'Reverse')
yticks([sqrt(2)-1, 0.5, 0.6, 0.7, 0.8, 0.9, 1])
yticklabels({'$\sqrt{2}-1$', '$0.5$', '$0.6$', '$0.7$', '$0.8$', '$0.9$', '$1$'})
xlabel('Normalized Acceleration')
ylabel('Normalized Escape $\Delta$v')
legend('Tangential', 'Azimuthal', 'location', 'northwest')

%print('DVesc.png', '-dpng', '-r600')

%% time

% figure('name', 't', 'units', 'inches', 'papersize', [6 6], 'paperposition', [0 0 6 6])
% grid on
% grid minor
% hold on
% plot(alphas, tesc_tan)
% plot(alphas, tesc_azi)
% ylim([10^-2, 10^5])
% xlim([alphas(end), alphas(1)])
% set(gca, 'Xscale', 'log')
% set(gca, 'Yscale', 'log')
% set(gca, 'Xdir', 'Reverse')
% xlabel('Normalized Acceleration')
% ylabel('Normalized Escape Burn Time')

%print('tesc.png', '-dpng', '-r600')

%% flight path angle

figure('name', 'gamma', 'units', 'inches', 'papersize', [6 6], 'paperposition', [0 0 6 6])
grid on
grid minor
hold on
plot(alphas, rad2deg(gammaesc_tan))
plot(alphas, rad2deg(gammaesc_azi))
ylim([0, 45])
xlim([alphas(end), alphas(1)])
set(gca, 'Xscale', 'log')
set(gca, 'Xdir', 'Reverse')
xlabel('Normalized Acceleration')
ylabel('Escape Flight Path Angle [deg]')
legend('Tangential', 'Azimuthal', 'location', 'northwest')

%print('gammaesc.png', '-dpng', '-r600')

%% delta v difference

figure
grid on
grid minor
hold on
plot(alphas, Dvesc_azi-Dvesc_tan, 'color', [0.4940, 0.1840, 0.5560])
set(gca, 'Xscale', 'log')
set(gca, 'Xdir', 'Reverse')
ylim([0 0.05])
yticks([0:0.01:0.05])
xlabel('Normalized Acceleration')
ylabel('Normalized Escape $\Delta$v Difference')

print('DvDelta.png', '-dpng', '-r600')

%%
