clc; clear; close all;

k = 1;
L = 1;
Nx = 100;
dx = L / Nx;
x = linspace(0, L, Nx+1);

dt = 0.00005;   % Stable time step
r = k * dt / dx^2;
fprintf("Fourier number r = %f\n", r);

if r > 0.5
    warning("FTCS is unstable because r > 0.5");
end

t_values = [0, 0.05, 0.10, 0.15];
Nt = max(round(t_values / dt));
t_steps = round(t_values / dt);

u = sin(pi * x);

U_plot = zeros(length(t_values), length(x));
U_plot(1, :) = u;

for n = 1:Nt
    u_new = u;
    
    for i = 2:Nx
        u_new(i) = u(i) + r * (u(i+1) - 2*u(i) + u(i-1));
    end
    
    u_new(1) = 0;
    u_new(end) = 0;
    u = u_new;

    for j = 2:length(t_values)
        if n == t_steps(j)
            U_plot(j, :) = u;
        end
    end
end

figure; hold on;
plot(x, U_plot(1,:), 'b', 'LineWidth', 2);
plot(x, U_plot(2,:), 'r', 'LineWidth', 2);
plot(x, U_plot(3,:), 'g', 'LineWidth', 2);
plot(x, U_plot(4,:), 'm', 'LineWidth', 2);

xlabel('x'); ylabel('u(x,t)');
title('Stable FTCS Numerical Solution');
legend('t=0','t=0.05','t=0.10','t=0.15');
grid on; hold off;
