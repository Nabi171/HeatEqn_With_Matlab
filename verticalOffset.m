clc; clear; close all;

% ============================
% Parameters
% ============================
k = 1;
L = 1;
Nx = 100;
dx = L / Nx;
x = linspace(0, L, Nx+1);

dt = 0.00005;
r = k * dt / dx^2;

t_compare = 0.15;
Nt = round(t_compare / dt);

% ============================
% Initial condition
% ============================
u_num = sin(pi * x);
u_exact = @(x,t) sin(pi*x).*exp(-k*pi^2*t);

% ============================
% FTCS Time-stepping
% ============================
for n = 1:Nt
    u_new = u_num;
    u_new(2:Nx) = u_num(2:Nx) + r*(u_num(3:Nx+1) - 2*u_num(2:Nx) + u_num(1:Nx-1));
    u_new(1) = 0;
    u_new(end) = 0;
    u_num = u_new;
end

% ============================
% Analytical & error
% ============================
u_ana = u_exact(x, t_compare);
rel_error = abs((u_ana - u_num) ./ u_ana);
rel_error(u_ana == 0) = 0;

% ============================
% Vertical Offsets
% ============================
offset1 = 1.2;      % analytical (top)
offset2 = 0.6;      % numerical (middle)
offset3 = 0.0;      % error (bottom)

% ============================
% Plot all curves in one space
% ============================
figure; hold on;

plot(x, u_ana + offset1, 'b', 'LineWidth', 2);        % Analytical shifted up
plot(x, u_num + offset2, 'r--', 'LineWidth', 2);      % Numerical shifted
plot(x, rel_error + offset3, 'k', 'LineWidth', 2);    % Relative error

title('Analytical, Numerical, and Relative Error (Separated Vertically)', 'FontSize', 16);
xlabel('x');
ylabel('Shifted Values');

% Custom y-ticks to show labels for each curve
yticks([offset1 offset2 0])
yticklabels({'Analytical', 'Numerical', 'Relative Error'})

grid on;
hold off;
