clc; clear; close all;

% ============================
% Parameters
% ============================
k = 1;                 % diffusion coefficient
L = 1;                 % domain length
Nx = 100;              % number of spatial steps
dx = L / Nx;           % spatial step size
x = linspace(0, L, Nx+1);

dt = 0.00005;          % stable time step (satisfies r <= 0.5)
r = k * dt / dx^2;
fprintf("\nFourier number r = %f\n", r);

T_final = 0.25;        % final simulation time
Nt = round(T_final / dt);

% ============================
% Initial and boundary conditions
% ============================
u = sin(pi * x);       % initial temperature
U_matrix = zeros(Nt+1, Nx+1);
U_matrix(1, :) = u;

% ============================
% FTCS Time Integration
% ============================
for n = 1:Nt
    u_new = u;
    for i = 2:Nx
        u_new(i) = u(i) + r * (u(i+1) - 2*u(i) + u(i-1));
    end

    % Boundary conditions
    u_new(1) = 0;
    u_new(end) = 0;

    u = u_new;
    U_matrix(n+1, :) = u;
end

% ============================
% Create 3D Surface Plot
% ============================
t = linspace(0, T_final, Nt+1);
[X, T] = meshgrid(x, t);

figure;
surf(X, T, U_matrix);

shading interp;
colormap turbo;
colorbar;

xlabel('Position (x)', 'FontSize', 14);
ylabel('Time (t)', 'FontSize', 14);
zlabel('Temperature u(x,t)', 'FontSize', 14);
title('3D Numerical FTCS Solution of the Heat Equation', 'FontSize', 16);

grid on;
view(45, 30);
