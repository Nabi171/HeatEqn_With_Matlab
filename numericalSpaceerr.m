clc; clear; close all;

% ============================
% Parameters
% ============================
k = 1;
L = 1;
Nx = 100;
dx = L / Nx;
x = linspace(0, L, Nx+1);

dt = 0.00005;                 % stable FTCS time step
r = k * dt / dx^2;

% Time at which to compare analytical & numerical
t_compare = 0.15;
Nt = round(t_compare / dt);

fprintf("Fourier number r = %f\n", r);

% ============================
% Initial condition
% ============================
u_num = sin(pi * x);          % numerical IC
u_exact = @(x,t) sin(pi*x).*exp(-k*pi^2*t);  % analytical function

% ============================
% FTCS Time-Stepping
% ============================
for n = 1:Nt
    u_new = u_num;

    u_new(2:Nx) = u_num(2:Nx) + r * (u_num(3:Nx+1) - 2*u_num(2:Nx) + u_num(1:Nx-1));

    u_new(1) = 0;
    u_new(end) = 0;

    u_num = u_new;
end

% ============================
% Analytical & Relative Error
% ============================
u_ana = u_exact(x, t_compare);
rel_error = abs((u_ana - u_num) ./ u_ana);
rel_error(u_ana == 0) = 0;

% ============================
% PLOTTING IN 3 SEPARATE SUBPLOTS
% ============================
figure;

% ---- Subplot 1: Analytical ----
subplot(3,1,1);
plot(x, u_ana, 'b', 'LineWidth', 2);
title('Analytical Solution at t = 0.15');
ylabel('u(x,t)');
grid on;

% ---- Subplot 2: Numerical ----
subplot(3,1,2);
plot(x, u_num, 'r--', 'LineWidth', 2);
title('Numerical FTCS Solution at t = 0.15');
ylabel('u(x,t)');
grid on;

% ---- Subplot 3: Relative Error ----
subplot(3,1,3);
plot(x, rel_error, 'k', 'LineWidth', 2);
title('Relative Error: |(Analytical - Numerical) / Analytical|');
xlabel('x');
ylabel('Error');
grid on;

sgtitle('Analytical, Numerical, and Relative Error Comparison');
