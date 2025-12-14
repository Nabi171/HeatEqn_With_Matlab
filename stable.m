clc; clear; close all;

k = 1;
L = 1;
Nx = 50;                   % Reduced for speed
dx = L / Nx;
x = linspace(0, L, Nx+1);

dt = 0.0001;               % Stable time step
r = k*dt/dx^2;
fprintf("Fourier number r = %f\n", r);

T_final = 0.15;            % Reduced final time (fast)
Nt = round(T_final/dt);

u = sin(pi*x);
U_matrix = zeros(Nt+1, Nx+1);
U_matrix(1,:) = u;

for n = 1:Nt
    u_new = u;
    u_new(2:Nx) = u(2:Nx) + r*(u(3:Nx+1) - 2*u(2:Nx) + u(1:Nx-1));

    u_new(1) = 0;
    u_new(end) = 0;

    u = u_new;
    U_matrix(n+1,:) = u;
end

t = linspace(0, T_final, Nt+1);
[X,T] = meshgrid(x,t);

figure;
surf(X,T,U_matrix);
shading interp;
colormap jet;
colorbar;

xlabel('x'); ylabel('t'); zlabel('u(x,t)');
title('3D Numerical FTCS Solution (Fast Version)');
view(45,30); grid on;
