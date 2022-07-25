function dSdt = odefun(t,S)
% The function was already rearranged by hand to solve for d(sigma)/dt
% Input equation with given variables

E_2A = 0.05;
K1 = 1.5;
K2 = 2.5;
eta = 100;
dE_1 = 0;

% For 2B:
% Given E = 0.1 * sin(t / 100), we can find derivative respect to t
dE_2 = 0.001 * cos(t / 100);
E_2B = 0.1 * sin(t/100);

% Set Function
dSdt(1) = ((dE_1 + ((K2/eta) * E_2A) - ((S(1) / eta)*((K1 + K2) / K1))) * K1);

dSdt = dSdt(:);

end
