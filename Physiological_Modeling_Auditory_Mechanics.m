% Michael Yang | my2699
% Due April 18, 2021
% BMENE6003
%% Homework 4B - Auditory Mechanics

clear; clc 
clear all;

%% Cochlea stimulated with pure-tone stimulation

% constants
p_0 = 2; %dyne / cm^2
% x_range = linspace(0, 3.5, 1000); %in cm
delx = 0.001;
x_range = [0:delx:3.5];

rho = 1; % g / cm^3
m = 0.01; % g/cm^2 - constant under x

% Stimulus Frequency
f = 2000; % Hz
omega = 2 * pi * f; % angular frequency

% Calculations
s = S(x_range);
r = R(x_range);
Z = 1i * ((m * omega) - (s / omega)) + r;
Y = 1./Z;
bracket = (-omega * rho * 2 * (cumsum(Y)*(delx)));
P_ans = p_0 * exp(bracket);

% Amplitude
P_amp = p_0 * exp(real(bracket));
% Phase
P_phase = imag(bracket);

% Question 1a: Initial model
figure(1);
subplot(1,3,1);
plot(x_range, P_amp);
title('BM Pressure Amplitude');
xlabel('cm');
ylabel('Pressure amp, dynes/cm^2');

subplot(1,3,2);
plot(x_range, (P_phase / (2*pi)));
title('BM Pressure Phase');
xlabel('cm');
ylabel('Pressure Phase, cycles');

subplot(1,3,3);
plot(x_range, P_ans);
title('BM Pressure Waveform');
xlabel('cm');
ylabel('pressure wave, dynes/cm^2');


% Question 1b: Velocity

%Velocity
BM_vel = 2 .* P_ans .* Y;

BM_vel_amp = abs(BM_vel);

BM_phase = unwrap(angle(BM_vel)) / (2*pi);

% Graphs for 1b
figure(2);
subplot(1,3,1);
plot(x_range, BM_vel_amp);
title('BM Velocity Amplitude');
xlabel('cm');
ylabel('BM Velocity amp, cm/s');

subplot(1,3,2);
plot(x_range, BM_phase);
title('BM Velocity Phase');
xlabel('cm');
ylabel('BM velocity phase, cycles');

subplot(1,3,3);
plot(x_range, BM_vel);
title('BM Velocity Waveform');
xlabel('cm');
ylabel('velocity wave, cm/s');


% 1c: Displacement

% Basic displacement
BM_disp = BM_vel / (omega);

% Displacement Waveform:
Disp_amp = abs(BM_disp);

% Displacement Phase
Disp_phase = unwrap(angle(BM_disp));

% Graphs for 1c
figure(3);
subplot(1,3,1);
plot(x_range, (Disp_amp * 1e7));
title('Displacement Amplitude');
xlabel('cm');
ylabel('Displacement amp, nm'); % MODIFY THIS

subplot(1,3,2);
plot(x_range, (Disp_phase / (2 * pi)));
title('Displacement Phase');
xlabel('cm');
ylabel('Disp Phase, cycles');

subplot(1,3,3);
plot(x_range, BM_disp);
title('Displacement Waveform');
xlabel('cm');
ylabel('Disp wave, cm');

% % 1d: Z
figure(4);
hold on;
plot(x_range, real(Z));
plot(x_range, imag(Z));
title('Real and Imaginary parts of Z');
xlabel('cm');
ylabel('dyne-s/cm^3');
legend('Real Z', 'Imaginary Z');

% 1e: Y
figure(5);
hold on;
plot(x_range, real(Y));
plot(x_range, imag(Y));
title('Real and Imaginary parts of Y');
xlabel('cm');
ylabel('1 / dyne-s/cm^3');
legend('Real Y', 'Imaginary Y');


%% Question 2
% constants
p_0 = 2; %dyne / cm^2
% x_range = linspace(0, 3.5, 1000); %in cm
delx = 0.001;
x_range = [0:delx:3.5];

rho = 1; % g / cm^3
m = 0.01; % g/cm^2 - constant under x

% Stimulus Frequency
f_range = [600 2000 6000]; % Hz
omega = 2 * pi .* f_range; % angular frequency

% Calculations

% 600 Hz
s = S(x_range);
r = R(x_range);
Z_600 = 1i * ((m * omega(1)) - (s / omega(1))) + r;
Y_600 = 1./Z_600;
bracket_600 = (-omega(1) * rho * 2 * (cumsum(Y_600)*(delx)));
P_ans_600 = p_0 * exp(bracket_600);

% 2000 Hz
Z_2000 = 1i * ((m * omega(2)) - (s / omega(2))) + r;
Y_2000 = 1./Z_2000;
bracket_2000 = (-omega(2) * rho * 2 * (cumsum(Y_2000)*(delx)));
P_ans_2000 = p_0 * exp(bracket_2000);

% 6000 Hz
Z_6000 = 1i * ((m * omega(3)) - (s / omega(3))) + r;
Y_6000 = 1./Z_6000;
bracket_6000 = (-omega(3) * rho * 2 * (cumsum(Y_6000)*(delx)));
P_ans_6000 = p_0 * exp(bracket_6000);

%Velocities: Amp + Phase
Stim_vel_600 = 2 .* P_ans_600 .* Y_600;
Stim_vel_2000 = 2 .* P_ans_2000 .* Y_2000;
Stim_vel_6000 = 2 .* P_ans_6000 .* Y_6000;

Stim_vel_amp_600 = abs(Stim_vel_600);
Stim_vel_amp_2000 = abs(Stim_vel_2000);
Stim_vel_amp_6000 = abs(Stim_vel_6000);

Stim_phase_600 = unwrap(angle(Stim_vel_600)) / (2*pi);
Stim_phase_2000 = unwrap(angle(Stim_vel_2000)) / (2*pi);
Stim_phase_6000 = unwrap(angle(Stim_vel_6000)) / (2*pi);


figure(6);
subplot(3,2,1);
plot(x_range, Stim_vel_amp_600);
title('V_{bm} (x) amplitude for 600 Hz');
xlabel('cm');
ylabel('cm/s');

subplot(3,2,2);
plot(x_range, Stim_phase_600);
title('V_{bm} (x) phase for 600 Hz');
xlabel('cm');
ylabel('cm/s');

subplot(3,2,3);
plot(x_range, Stim_vel_amp_2000);
title('V_{bm} (x) amplitude for 2000 Hz');
xlabel('cm');
ylabel('cm/s');

subplot(3,2,4);
plot(x_range, Stim_phase_2000);
title('V_{bm} (x) phase for 2000 Hz');
xlabel('cm');
ylabel('cm/s');

subplot(3,2,5);
plot(x_range, Stim_vel_amp_6000);
title('V_{bm} (x) amplitude for 6000 Hz');
xlabel('cm');
ylabel('cm/s');

subplot(3,2,6);
plot(x_range, Stim_phase_6000);
title('V_{bm} (x) phase for 6000 Hz');
xlabel('cm');
ylabel('cm/s');


%% Question 3

% COND 1: m = 0
% constants
p_0 = 2; %dyne / cm^2
% x_range = linspace(0, 3.5, 1000); %in cm
delx = 0.001;
x_range = [0:delx:3.5];

rho = 1; % g / cm^3
m = 0; % g/cm^2 - CHANGED

% Stimulus Frequency
f = 2000; % Hz
omega = 2 * pi * f; % angular frequency

% Calculations
s = S(x_range);
r = R(x_range);
Z_m = 1i * ((m * omega) - (s / omega)) + r;
Y_m = 1./Z_m;
bracket_m = (-omega * rho * 2 * (cumsum(Y_m)*(delx)));
P_ans_m = p_0 * exp(bracket_m);

% Velocity (M Change): Amplitude and Phase
Vel_m = 2 .* P_ans_m .* Y_m;

m_vel_amp = abs(Vel_m);

m_phase = unwrap(angle(Vel_m)) / (2*pi);

% Graph change with m
figure(7);
subplot(1,2,1);
plot(x_range, m_vel_amp);
title('V_{bm} (x) amplitude for m = 0');
xlabel('cm');
ylabel('cm/s');

subplot(1,2,2);
plot(x_range, m_phase);
title('V_{bm} (x) phase for m = 0');
xlabel('cm');
ylabel('cm/s');
%------------------------------------------------------

% COND 2: r * 0.1
% constants
p_0 = 2; %dyne / cm^2
% x_range = linspace(0, 3.5, 1000); %in cm
delx = 0.001;
x_range = [0:delx:3.5];

rho = 1; % g / cm^3
m = 0.01; % g/cm^2 - BACK TO NORMAL

% Stimulus Frequency
f = 2000; % Hz
omega = 2 * pi * f; % angular frequency

% Calculations
s = S(x_range);
r = R(x_range) .* 0.1; % - CHANGED
Z_r = 1i * ((m * omega) - (s / omega)) + r;
Y_r = 1./Z_r;
bracket_r = (-omega * rho * 2 * (cumsum(Y_r)*(delx)));
P_ans_r = p_0 * exp(bracket_r);

% Velocity (M Change): Amplitude and Phase
Vel_r = 2 .* P_ans_r .* Y_r;

r_vel_amp = abs(Vel_r);

r_phase = unwrap(angle(Vel_r)) / (2*pi);

% Graph change with R
figure(8);
subplot(1,2,1);
plot(x_range, r_vel_amp);
title('V_{bm} (x) amplitude for 10% of R');
xlabel('cm');
ylabel('cm/s');

subplot(1,2,2);
plot(x_range, r_phase);
title('V_{bm} (x) phase for 10% of R');
xlabel('cm');
ylabel('cm/s');
%------------------------------------------------------

% COND 3: s * 10
% constants
p_0 = 2; %dyne / cm^2
% x_range = linspace(0, 3.5, 1000); %in cm
delx = 0.001;
x_range = [0:delx:3.5];

rho = 1; % g / cm^3
m = 0.01; % g/cm^2 - BACK TO NORMAL

% Stimulus Frequency
f = 2000; % Hz
omega = 2 * pi * f; % angular frequency

% Calculations
s = S(x_range) .* 10; % - CHANGED
r = R(x_range); 
Z_s = 1i * ((m * omega) - (s / omega)) + r;
Y_s = 1./Z_s;
bracket_s = (-omega * rho * 2 * (cumsum(Y_s)*(delx)));
P_ans_s = p_0 * exp(bracket_s);

% Velocity (M Change): Amplitude and Phase
Vel_s = 2 .* P_ans_s .* Y_s;

s_vel_amp = abs(Vel_s);

s_phase = unwrap(angle(Vel_s)) / (2*pi);

% Graph change with S
figure(9);
subplot(1,2,1);
plot(x_range, s_vel_amp);
title('V_{bm} (x) amplitude for 10 x S');
xlabel('cm');
ylabel('cm/s');

subplot(1,2,2);
plot(x_range, s_phase);
title('V_{bm} (x) phase for 10 x S');
xlabel('cm');
ylabel('cm/s');
%-------------------------------------------------------

%% Question 4

% set point 1.5 cm
pt = 1.5;

% omega
omega = 2 * pi * 2000; % Hz

% Finding m at s/w^2
S_pt = S(pt);
M_pt = S_pt / (omega^2);

% older constants
p_0 = 2; %dyne / cm^2
% x_range = linspace(0, 3.5, 1000); %in cm
delx = 0.001;
x_range = [0:delx:3.5];

rho = 1; % g / cm^3

% Calculations - make S normal here
s = S(x_range);
r = R(x_range); 
Z_4 = 1i * ((M_pt * omega) - (s / omega)) + r;
Y_4 = 1./Z_4;
bracket_4 = (-omega * rho * 2 * (cumsum(Y_4)*(delx)));
P_ans_4 = p_0 * exp(bracket_4);

% Velocity (M Change): Amplitude and Phase
Vel_4 = 2 .* P_ans_4 .* Y_4;

vel_amp_4 = abs(Vel_4);

phase_4 = unwrap(angle(Vel_4)) / (2*pi);

% Plot
figure(10);
subplot(1,2,1);
plot(x_range, vel_amp_4);
title('V_{bm} (x) amplitude with M at local resonance');
xlabel('cm');
ylabel('cm/s');

subplot(1,2,2);
plot(x_range, phase_4);
title('V_{bm} (x) phase with M at local resonance');
xlabel('cm');
ylabel('cm/s');

% Z 
figure(11);
hold on;
plot(x_range, real(Z_4));
plot(x_range, imag(Z_4));
title('Real and Imaginary parts of Z, R normal');
xlabel('cm');
ylabel('dyne-s/cm^3');
legend('Real Z', 'Imaginary Z');

%--------------------------------------------------------
% Q4c - R * 0.1
% Calculations
r = R(x_range) * 0.1; 
Z_4c = 1i * ((M_pt * omega) - (s / omega)) + r;
Y_4c = 1./Z_4c;
bracket_4c = (-omega * rho * 2 * (cumsum(Y_4c)*(delx)));
P_ans_4c = p_0 * exp(bracket_4c);

% Velocity (M Change): Amplitude and Phase
Vel_4c = 2 .* P_ans_4c .* Y_4c;

vel_amp_4c = abs(Vel_4c);

phase_4c = unwrap(angle(Vel_4c)) / (2*pi);

% Plot
figure(12);
subplot(1,2,1);
plot(x_range, vel_amp_4c);
title('V_{bm} (x) amplitude with M changed, R changed');
xlabel('cm');
ylabel('cm/s');

subplot(1,2,2);
plot(x_range, phase_4c);
title('V_{bm} (x) phase with M changed, R changed');
xlabel('cm');
ylabel('cm/s');

% Z 
figure(13);
hold on;
plot(x_range, real(Z_4c));
plot(x_range, imag(Z_4c));
title('Real and Imaginary parts of Z, R changed');
xlabel('cm');
ylabel('dyne-s/cm^3');
legend('Real Z', 'Imaginary Z');

%--------------------------------------------------

%% Q6: Changing Impedance by Adding Negative Resistance

% omega
omega = 2 * pi * 2000; % Hz

% older constants
p_0 = 2; %dyne / cm^2
delx = 0.001;
x_range = [0:delx:3.5];

rho = 1; % g / cm^3
m = 0.01;

% Establish R so it can be selectively changed
r_amp = R(x_range);

% Calculations: changing R in range
for i = 1:length(x_range)
    s = S(x_range);
    if (1 < x_range(i)) && (x_range(i) < 2)
        % Introduce negative resistance for short x
        r_amp(i) = -5*exp(2.25 * x_range(i));
    end
end

% Rest of calculations
Z_5 = 1i * ((m * omega) - (s / omega)) + r_amp;
Y_5 = 1./Z_5;
bracket_5 = (-omega * rho * 2 * (cumsum(Y_5)*(delx)));
P_ans_5 = p_0 * exp(bracket_5);

% Velocity with Sitmulus Enhancement
Vel_5 = 2 .* P_ans_5 .* Y_5;

vel_amp_5 = abs(Vel_5);

V_phase_5 = unwrap(angle(Vel_5)) / (2*pi);

% Pressure Amp + Phase
% Amplitude
P_amp_5 = p_0 * exp(real(bracket_5));
% Phase
P_phase_5 = imag(bracket_5);

figure(14);
hold on;
plot(x_range, real(Z_5));
plot(x_range, imag(Z_5));
title('Real and Imaginary parts of Z, Amplifier');
xlabel('cm');
ylabel('dyne-s/cm^3');
legend('Real Z', 'Imaginary Z');

% Velocity Amp + Phase
figure(15);
subplot(1,2,1);
plot(x_range, vel_amp_5);
title('V_{bm} (x) amplitude with Stimulus Enhancement');
xlabel('cm');
ylabel('cm/s');

subplot(1,2,2);
plot(x_range, V_phase_5);
title('V_{bm} (x) phase with Stimulus Enhancement');
xlabel('cm');
ylabel('cm/s');


% Pressure Amp + Phase
figure(16);
subplot(1,2,1);
plot(x_range, P_amp_5);
title('P(x) Amplitude with Stimulus Enhancement');
xlabel('cm');
ylabel('Pressure amp, dynes/cm^2');

subplot(1,2,2);
plot(x_range, (P_phase_5 / (2*pi)));
title('P(x) Phase with Stimulus Enhancement');
xlabel('cm');
ylabel('Pressure Phase, cycles');

%% Question 6: Stiffness or Damping

% omega
omega = 2 * pi * 2000; % Hz

% older constants
p_0 = 2; %dyne / cm^2
delx = 0.001;
x_range = [0:delx:3.5];

rho = 1; % g / cm^3
m = 0.01;

% Establish S so it can be selectively changed
s_change = S(x_range);

% Calculations: changing S in beginning range
for i = 1:length(x_range)
    r = R(x_range);
    if (0 < x_range(i)) && (x_range(i) < 1)
        % Introduce negative resistance for short x
        s_change(i) = 2e8*exp(1.5 * x_range(i));
        % Or, for no dependence
%         s_change(i) = 5e8;
    end
end

% Rest of calculations
Z_6 = 1i * ((m * omega) - (s_change / omega)) + r;
Y_6 = 1./Z_6;
bracket_6 = (-omega * rho * 2 * (cumsum(Y_6)*(delx)));
P_ans_6 = p_0 * exp(bracket_6);

% Velocity with Sitmulus Enhancement
Vel_6 = 2 .* P_ans_6 .* Y_6;

vel_amp_6 = abs(Vel_6);

V_phase_6 = unwrap(angle(Vel_6)) / (2*pi);

% Pressure Amp + Phase
% Amplitude
P_amp_6 = p_0 * exp(real(bracket_6));
% Phase
P_phase_6 = imag(bracket_6);

figure(17);
hold on;
plot(x_range, real(Z_6));
plot(x_range, imag(Z_6));
title('Real and Imaginary parts of Z, Amplifier');
xlabel('cm');
ylabel('dyne-s/cm^3');
legend('Real Z', 'Imaginary Z');

% Velocity Amp + Phase
figure(18);
subplot(1,2,1);
plot(x_range, vel_amp_6);
title('V_{bm} (x) amplitude with Implant');
xlabel('cm');
ylabel('cm/s');

subplot(1,2,2);
plot(x_range, V_phase_6);
title('V_{bm} (x) phase with Implant');
xlabel('cm');
ylabel('cm/s');


% Pressure Amp + Phase
figure(19);
subplot(1,2,1);
plot(x_range, P_amp_6);
title('P(x) Amplitude with Implant');
xlabel('cm');
ylabel('Pressure amp, dynes/cm^2');

subplot(1,2,2);
plot(x_range, (P_phase_6 / (2*pi)));
title('P(x) Phase with Stimulus Enhancement');
xlabel('cm');
ylabel('Pressure Phase, cycles');



% Functions
function calc_s = S(x)
    calc_s = 2e8 * exp(-1.5 * x); % returns in dyne/cm^3
end

function calc_r = R(x)
    calc_r = 5*exp(2.25 * x); % dynes-s/cm^3
end

% function calc_p = P(p_0, omega, rho, x, y)
%     calc_p = p_0 * exp(-2 * omega * rho * cumsum(y)*x);
% end


% P = (P_0) * e^(1i * w * t) * e^(-1i * k * x_range) * e^(-k * z);