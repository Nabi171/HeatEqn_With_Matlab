% 3D Heat equation - analytical series solution
% u_t = k(u_xx + u_yy + u_zz)

clear; close all; clc;

% Parameters
Lx = 1; Ly = 1; Lz = 1;
k  = 1;
M = 10; N = 10; P = 10;   % number of modes

% Spatial grid
Nx = 60; Ny = 60; Nz = 60;
x = linspace(0,Lx,Nx);
y = linspace(0,Ly,Ny);
z = linspace(0,Lz,Nz);

[X,Y,Z] = meshgrid(x,y,z);

% Time
t = 0.05;

% Initial condition (single eigenmode example)
u0 = @(x,y,z) sin(pi*x/Lx).*sin(pi*y/Ly).*sin(pi*z/Lz);

% Compute coefficients analytically (only one nonzero here)
A = zeros(M,N,P);
A(1,1,1) = 1;  % because u0 is exactly first eigenmode

% Build solution
U = zeros(size(X));

for l = 1:M
    for m = 1:N
        for n = 1:P
            lambda = k*pi^2*(l^2/Lx^2 + m^2/Ly^2 + n^2/Lz^2);
            U = U + A(l,m,n) .* ...
                sin(l*pi*X/Lx) .* ...
                sin(m*pi*Y/Ly) .* ...
                sin(n*pi*Z/Lz) .* ...
                exp(-lambda*t);
        end
    end
end

% Plot a slice at z = Lz/2
z_index = round(Nz/2);
figure;
surf(X(:,:,z_index), Y(:,:,z_index), U(:,:,z_index), 'EdgeColor','none');
xlabel('x'); ylabel('y'); zlabel('u');
title('3D Heat Equation Solution (z = L_z/2 slice)');
colorbar;
view(35,30);
