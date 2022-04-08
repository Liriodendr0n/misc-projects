clear
close all
clc

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

load foils.mat


sparXs = [sin(pi/10)^2 sin(3*pi/10)^2];


% %% 13%
% 
% coords = [af13.x.*sign(af13.y), af13.y];
% coords = coords(2:end, :);
% 
% for i = 1:length(sparXs)
%     sparYs(i,:) = [interp1(coords(:,1), coords(:,2), sparXs(i)); ...
%         interp1(coords(:,1), coords(:,2), -sparXs(i))];
%     
% end
% 
% 
% figure
% hold on
% axis equal
% plot(af13.x, af13.y)
% % plot(coords(:,1), coords(:,2))
% for i = 1:length(sparXs)
%     plot([sparXs(i), sparXs(i)], sparYs(i,:))
% end
% 
% sparHs(1,:) = -diff(sparYs,1,2);

%% 17%

coords = [af17.x.*sign(af17.y), af17.y];
coords = coords(2:end, :);

for i = 1:length(sparXs)
    sparYs(i,:) = [interp1(coords(:,1), coords(:,2), sparXs(i)); ...
        interp1(coords(:,1), coords(:,2), -sparXs(i))];
    
end


figure('units', 'inches', 'papersize', [7 3], 'paperposition', [0 0 7 3])
hold on
grid on
grid minor
axis equal
plot(af17.x, af17.y, 'color', 'black')
% plot(coords(:,1), coords(:,2))
for i = 1:length(sparXs)
    plot([sparXs(i), sparXs(i)], sparYs(i,:), 'linewidth', 2)
end

sparHs(2,:) = -diff(sparYs,1,2)
xlim([-0.1 1.1])
ylim([-0.2 0.2])
xlabel('x/c')
ylabel('y/c')
xticks([0 0.15 0.3 0.5 0.65 1])

print('sparLocs.pdf', '-dpdf', '-painters', '-fillpage')

%%




