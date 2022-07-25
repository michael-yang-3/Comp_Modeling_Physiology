% Michael Yang | my2699
% Due April 4, 2021
% BMENE6003
%% Homework 4A - Auditory Mechanics

clear; clc 
clear all;

%% Question 7A-B

% Impedance Z = i * (m*w - s/w) + r

% Constants
m = 0.015; % g/cm^2
s = 2*(10^7); % dyne/cm^3, or g/(s^2 cm^2)
r = 600; %g/s*cm^2

% Check based on class
m2 = 0.1;
s2 = 5 * (10^7);
r2 = 400;

% Plot along frequencies of 100 Hz to 10 kHz
W_range = (linspace(100,10000,100));

%Frequency = w / 2*pi --> w = frequency * 2*pi

Z_ans = zeros(length(W_range));

for i = 1 : length(W_range)
    w_step = W_range(i) * 2 * pi;
    x = Calc_Imp(m, w_step, s, r);
    Z_ans(i) = x;
end

Z_ans_real = real(Z_ans);
Z_ans_imag = imag(Z_ans);

% Question 7a
win = figure(1);
win(1) = subplot(2,3,1);
set(win, 'Nextplot', 'add');
plot(win(1), W_range, Z_ans_imag);
plot(win(1), W_range, Z_ans_real);

xlim([100 10000]);
ylim([-inf 2e4]);
title('Impedance versus Frequency');
xlabel('Frequency (Hz)');
ylabel('Impedance');

legend('Imaginary Impedance', 'Real Impedance');

% Question 7B
% Magnitude
Z_mag = zeros(length(W_range));
Z_theta = zeros(length(W_range));

for i = 1 : length(W_range)
    w_step = W_range(i) * 2 * pi;
    x = Calc_Imp(m, w_step, s, r);
    Z_ans(i) = x;
    Z_mag(i) = (((real(x))^2 + (imag(x))^2)^0.5);
    Z_theta(i) = (atan(imag(x) / real(x)));
end

subplot(2,3,2);
hold on;
plot(W_range, Z_mag);
title('Impedance Mag versus Frequency');
xlabel('Frequency (Hz)');
ylabel('Impedance Magnitude (g/(s*cm^2))');

% Phase
subplot(2,3,3);
hold on;
plot(W_range, Z_theta);
title('Impedance Phase versus Frequency');
xlabel('Frequency (Hz)');
ylabel('Z Phase, in Radians');


% Question 7C-D
Y_admit = zeros(length(W_range));

Y_mag = zeros(length(W_range));
Y_theta = zeros(length(W_range));

% Loop over frequencies
for i = 1 : length(W_range)
    w_step = W_range(i) * 2 * pi;
    x = Calc_Imp(m, w_step, s, r);
    Z_ans(i) = x;
    Y_admit(i) = (1 / x);
    Y_mag(i) = (((real(Y_admit(i)))^2 + (imag(Y_admit(i)))^2)^0.5);
    Y_theta(i) = (atan(imag(Y_admit(i)) / real(Y_admit(i))));
end

Y_ans_real = real(Y_admit);
Y_ans_imag = imag(Y_admit);

% Admittance (real and imaginary) vs frequency
subplot(2,3,4);
hold on;
plot(W_range, Y_ans_real); 
plot(W_range, Y_ans_imag);

xlim([100 10000]);
ylim([-inf 2.5e-3]);
title('Admittance versus Frequency');
xlabel('Frequency (Hz)');
ylabel('Admittance (real and imaginary)');

legend('Real Admittance', 'Imaginary Admittance');

% Question 7D - Admittance Magnitude and Phase
subplot(2,3,5);
hold on; 
plot(W_range, Y_mag);
title('Admittance Mag versus Frequency');
xlabel('Frequency (Hz)');
ylabel('Admittance Magnitude (s*cm^2)/g');

% Phase
subplot(2,3,6);
hold on; 
plot(W_range, Y_theta);
title('Admittance Phase versus Frequency');
xlabel('Frequency (Hz)');
ylabel('Y Phase, in Radians');

% Function
function imp = Calc_Imp(m, w, s, r)
    imp = 1i * ((m*w) - (s / w)) + r;
end
