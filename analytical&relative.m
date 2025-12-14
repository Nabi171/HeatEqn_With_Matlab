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
% Analytical + Error
% ============================
u_ana = u_exact(x, t_compare);
rel_error = abs((u_ana - u_num) ./ u_ana);
rel_error(u_ana == 0) = 0;

% ============================
% OFFSETS
% ============================
offset_error = -0.8;   % Put error curve lower

% ============================
% PLOT
% ============================
figure; hold on;

% Analytical and Numerical almost overlapping
plot(x, u_ana, 'b', 'LineWidth', 2);
plot(x, u_num, 'r--', 'LineWidth', 2);

% Error curve shifted downward
plot(x, rel_error + offset_error, 'k', 'LineWidth', 2);

title('Analytical, Numerical, and Relative Error (Close + Offset Layout)');
xlabel('x');
ylabel('Values');

% Label y-axis ticks
yticks([0, offset_error])
yticklabels({'Analytical & Numerical', 'Relative Error'});

legend('Analytical','Numerical','Relative Error (shifted)', ...
       'Location','northeastoutside');

grid on;
hold off;
