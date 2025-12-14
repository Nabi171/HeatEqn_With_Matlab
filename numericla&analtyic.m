clc; clear; close all;

% ============================
% Parameters
% ============================
k = 1;
L = 1;
Nx = 100;                % spatial steps
dx = L / Nx;
x = linspace(0, L, Nx+1);

dt = 0.00005;            % stable time step (r <= 0.5)
r = k * dt / dx^2;

% Select a time for comparison
t_compare = 0.15;        % time level for error analysis
Nt = round(t_compare / dt);

fprintf("Fourier number r = %f\n", r);

% ============================
% Initial Condition
% ============================
u_num = sin(pi * x);     % numerical initial condition
u_exact = @(x,t) sin(pi*x).*exp(-k*pi^2*t);

% ============================
% FTCS Time-stepping
% ============================
for n = 1:Nt
    u_new = u_num;
    
    % Update interior points using FTCS
    u_new(2:Nx) = u_num(2:Nx) + r * (u_num(3:Nx+1) - 2*u_num(2:Nx) + u_num(1:Nx-1));
    
    % Boundary conditions
    u_new(1) = 0;
    u_new(end) = 0;
    
    u_num = u_new;
end

% ============================
% Analytical Solution at t = t_compare
% ============================
u_ana = u_exact(x, t_compare);

% ============================
% Relative Error Calculation
% ============================
rel_error = abs((u_ana - u_num) ./ u_ana);
rel_error(u_ana == 0) = 0; % prevent NaN near boundaries

% ============================
% Combined Plot
% ============================
figure;
hold on;

plot(x, u_ana, 'b', 'LineWidth', 2);
plot(x, u_num, 'r--', 'LineWidth', 2);
plot(x, rel_error, 'k', 'LineWidth', 2);

xlabel('x', 'FontSize', 14);
ylabel('Value', 'FontSize', 14);
title('Analytical, Numerical, and Relative Error Plot of u(x,t)', 'FontSize', 16);

legend('Analytical Solution', 'Numerical FTCS', 'Relative Error', ...
    'Location','northeastoutside');

grid on;
hold off;
