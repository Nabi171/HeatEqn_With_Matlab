clc; clear; close all;

% ============================
% Parameters
% ============================
k = 1;                 % Thermal diffusivity
L = 1;                 % Domain length (0 to 1)
Nx = 100;              % Number of spatial divisions
dx = L / Nx;           % Spatial step size
x = linspace(0, L, Nx+1);

dt = 0.00005;          % Stable time step (FTCS requires r <= 0.5)
r = k * dt / dx^2;     % Fourier number

fprintf("\nFourier number r = %f\n", r);

% ============================
% Time configuration
% ============================
T_final = 0.25;               % Final simulation time
Nt = round(T_final / dt);     % Number of time steps

% ============================
% Initial condition u(x,0)
% ============================
u = sin(pi * x);              % Initial temperature profile

% ============================
% Storage matrix for 3D plotting
% ============================
U_matrix = zeros(Nt+1, Nx+1);
U_matrix(1, :) = u;           % Store t=0

% ============================
% FTCS Time-stepping Loop
% ============================
for n = 1:Nt
    u_new = u;

    % Update interior nodes using FTCS formula
    for i = 2:Nx
        u_new(i) = u(i) + r * (u(i+1) - 2*u(i) + u(i-1));
    end

    % Boundary Conditions: u(0,t) = 0 and u(1,t) = 0
    u_new(1) = 0;
    u_new(end) = 0;

    u = u_new;
    U_matrix(n+1, :) = u;     % Store solution at this time step
end

% ============================
% Prepare mesh for 3D surface
% ============================
t = linspace(0, T_final, Nt+1);
[X, T] = meshgrid(x, t);

% ============================
% Plot 3D Numerical FTCS Solution
% ============================
figure;
surf(X, T, U_matrix);

shading interp;               % Smooth shading
colormap jet;                 % Jet colormap for temperature variation
colorbar;                     % Add color legend

xlabel('Position (x)', 'FontSize', 14);
ylabel('Time (t)', 'FontSize', 14);
zlabel('Temperature u(x,t)', 'FontSize', 14);
title('3D Numerical FTCS Solution of the Heat Equation', 'FontSize', 16);

grid on;
view(45, 30);                 % 3D viewing angle
