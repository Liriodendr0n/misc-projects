clear
close all
clc

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% description
% this is a brute force approach to the sizing problem
% I choose a parameter space, and try every member of that space
% configurations that meet every constraint get scored
% the configuration with the highest score wins
% kriging or co-kriging are the right ways to do this, I will not use them

%% define parameter space
tic
parameters.S = linspace(2, 16, 15);
parameters.b = linspace(2, 8, 13);
parameters.Nbox = linspace(1, 20, 20);
parameters.CLmax = 1.5;
parameters.EWF = 0.7;
parameters.P = 1200;

% this will be very large
Ncases = prod(structfun(@numel, parameters));

%% permute cases

% form multidimensional array of ones of the correct shape
tempidx = ones(structfun(@numel, parameters)');
% get indices of every element in the array
[idx.S, idx.b, idx.Nbox, idx.CLmax, idx.EWF, idx.P] = findND(tempidx);

% permute every configuration
for i = 1:Ncases
    configuration(i).S = parameters.S(idx.S(i));
    configuration(i).b = parameters.b(idx.b(i));
    configuration(i).Nbox = parameters.Nbox(idx.Nbox(i));
    configuration(i).CLmax = parameters.CLmax(idx.CLmax(i));
    configuration(i).EWF = parameters.EWF(idx.EWF(i));
    configuration(i).P = parameters.P(idx.P(i));
end

%% evaluate configurations


[configuration] = evalConfig(configuration);

for i = 1:Ncases
    isValid(i) = configuration(i).isValid;
end

valids = find(isValid);

[M2val, M2id] = max([configuration(valids).M2score]);
[M3val, M3id] = max([configuration(valids).M3score]);

%% display winner
disp(" ")
disp("Success Fraction: " + length(valids)/Ncases)
disp("Winning Configuration: " + valids(M2id))
disp(configuration(valids(M2id)))

%% plots

toc
% % wing area
% N = 301;
% xv = linspace(min([configuration(valids).M2score]), max([configuration(valids).M2score]), N);
% yv = linspace(min([configuration(valids).M3score]), max([configuration(valids).M3score]), N);
% 
% [X, Y] = meshgrid(xv, yv);
% Z = griddata([configuration(valids).M2score], [configuration(valids).M3score], [configuration(valids).S], X, Y);
% 
% figure('name', 'Strend', 'papersize', [6 4])
% hold on
% grid on
% contourf(X, Y, Z, 30, 'linecolor', 'none')
% title('Wing Area Trends')
% xlabel('Mission 2 Score')
% ylabel('Mission 3 Score')
% xlim([0 4])
% ylim([0 20])
% c = colorbar;
% c.Label.String = '$S$';
% c.Label.Interpreter = 'latex';
% %print('Strend.png', '-dpng', '-r600')
% 
% 
% 
% % wing span
% N = 301;
% xv = linspace(min([configuration(valids).M2score]), max([configuration(valids).M2score]), N);
% yv = linspace(min([configuration(valids).M3score]), max([configuration(valids).M3score]), N);
% 
% [X, Y] = meshgrid(xv, yv);
% Z = griddata([configuration(valids).M2score], [configuration(valids).M3score], [configuration(valids).b], X, Y);
% 
% figure('name', 'btrend', 'papersize', [6 4])
% hold on
% grid on
% contourf(X, Y, Z, 30, 'linecolor', 'none')
% title('Wing Span Trends')
% xlabel('Mission 2 Score')
% ylabel('Mission 3 Score')
% xlim([0 4])
% ylim([0 20])
% c = colorbar;
% c.Label.String = '$b$';
% c.Label.Interpreter = 'latex';
% %print('btrend.png', '-dpng', '-r600')
% 
% 
% % CLmax
% N = 301;
% xv = linspace(min([configuration(valids).M2score]), max([configuration(valids).M2score]), N);
% yv = linspace(min([configuration(valids).M3score]), max([configuration(valids).M3score]), N);
% 
% [X, Y] = meshgrid(xv, yv);
% Z = griddata([configuration(valids).M2score], [configuration(valids).M3score], [configuration(valids).CLmax], X, Y);
% 
% figure('name', 'CLmaxtrend', 'papersize', [6 4])
% hold on
% grid on
% contourf(X, Y, Z, 30, 'linecolor', 'none')
% title('CLmax Trends')
% xlabel('Mission 2 Score')
% ylabel('Mission 3 Score')
% xlim([0 4])
% ylim([0 20])
% c = colorbar;
% c.Label.String = '$C_{Lmax}$';
% c.Label.Interpreter = 'latex';
% %print('CLmaxtrend.png', '-dpng', '-r600')
% 
% 
% 
% % P
% N = 301;
% xv = linspace(min([configuration(valids).M2score]), max([configuration(valids).M2score]), N);
% yv = linspace(min([configuration(valids).M3score]), max([configuration(valids).M3score]), N);
% 
% [X, Y] = meshgrid(xv, yv);
% Z = griddata([configuration(valids).M2score], [configuration(valids).M3score], [configuration(valids).P], X, Y);
% 
% figure('name', 'Ptrend', 'papersize', [6 4])
% hold on
% grid on
% contourf(X, Y, Z, 30, 'linecolor', 'none')
% title('Power Trends')
% xlabel('Mission 2 Score')
% ylabel('Mission 3 Score')
% xlim([0 4])
% ylim([0 20])
% c = colorbar;
% c.Label.String = '$P$';
% c.Label.Interpreter = 'latex';
% %print('TWtrend.png', '-dpng', '-r600')
% 



