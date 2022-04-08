clear
clc
close all

load wingfoils.mat

figure
subplot(2,1,1)
hold on
plot(MS10313.x, MS10313.y)
plot(MSc0613.x, MSc0613.y)
axis equal

subplot(2,1,2)
hold on
plot(MS10317.x, MS10317.y)
plot(MSc0617.x, MSc0617.y)
axis equal



