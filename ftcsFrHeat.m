clc; clear; close all;

% ============================
% Parameters
% ============================
k = 1;                 % Thermal diffusivity
L = 1;                 % Length of the rod
Nx = 100;              % Number of spatial divisions
dx = L / Nx;           % Spatial step size
x = linspace(0, L, Nx+1);

dt = 0.0001;           % Time step (choose small enough for stability)
r = k * dt / dx^2;     % Fourier number

fprintf("Fourier number r = %f\n", r);

% Must satisfy FTCS stability: r <= 0.5
if r > 0.5
    warning("FTCS is unstable because r > 0.5");
end

% ============================
% Time selection for plotting
% ============================
t_values = [0, 0.05, 0.10, 0.15];   % times of interest
Nt = max(round(t_values / dt));     % total steps needed
t_steps = round(t_values / dt);     % convert to index steps

% ============================
% Initial condition
% ============================
u = sin(pi * x);       % u(x,0)
U_plot = zeros(length(t_values), length(x));
U_plot(1, :) = u;      % Store t = 0

% ============================
% FTCS Time-stepping
% ============================
for n = 1:Nt
    u_new = u;
    
    % Update interior nodes
    for i = 2:Nx
        u_new(i) = u(i) + r * (u(i+1) - 2*u(i) + u(i-1));
    end
    
    % Boundary conditions
    u_new(1) = 0;
    u_new(end) = 0;

    % Update the solution
    u = u_new;
    
    % Save selected time snapshots
    for j = 2:length(t_values)
        if n == t_steps(j)
            U_plot(j, :) = u;
        end
    end
end

% ============================
% Plot Numerical FTCS Solution
% ============================
figure;
hold on;

plot(x, U_plot(1,:), 'b', 'LineWidth', 2);
plot(x, U_plot(2,:), 'r', 'LineWidth', 2);
plot(x, U_plot(3,:), 'g', 'LineWidth', 2);
plot(x, U_plot(4,:), 'm', 'LineWidth', 2);

xlabel('x', 'FontSize', 14);
ylabel('u(x,t)', 'FontSize', 14);
title('Numerical FTCS Solution of the Heat Equation', 'FontSize', 16);

legend('t = 0', 't = 0.05', 't = 0.10', 't = 0.15', ...
       'Location', 'northeastoutside');

grid on;
hold off;
