% Michael Yang | my2699
% Due Feb 14, 2022
% BMEN E6003
%% Homework 2A Problem 2

% to run your code simply hit "Run" or use the shortcut "command(Mac)/
% control(Windows)+enter"

clear;clc

% Set t-span, in seconds
tspan = 1800;

% Initial Conditions
S_0A = 0.25;
% S_0B = 0;
y0 = [S_0A];
% y0 = [S_0B];

% run ode23t
[time, S] = ode23t(@odefun, [0 tspan], y0);

% Graph for t-span of the overall stress over time
figure;
subplot(1, 1, 1);
plot(time, S(:));
title('Overall Stress over Time - Stress-Relaxation');
xlabel('Time (seconds)');
ylabel('Stress (MPa)');

% plot 2B
% figure;
% subplot(1, 1, 1);
% plot(time, S(:));
% title('Overall Stress over Time - Sinusoidal');
% xlabel('Time (seconds)');
% ylabel('Stress (MPa)');