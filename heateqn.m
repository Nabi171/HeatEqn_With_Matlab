% 2D Heat equation - analytical series solution (Dirichlet 0 on boundary)
% u_t = k (u_xx + u_yy), domain (0,Lx)x(0,Ly)
clear; close all; clc;

% Parameters
Lx = 1; Ly = 1;      % domain sizes
k  = 1;              % thermal diffusivity
M  = 30; N = 30;     % number of eigenmodes in x and y (increase for accuracy)

% Spatial grid for evaluating & plotting
Nx = 80; Ny = 80;
x = linspace(0,Lx,Nx);
y = linspace(0,Ly,Ny);
[X,Y] = meshgrid(x,y);

% Time points to plot
t_plot = [0, 0.01, 0.05, 0.2];

% Initial condition u0(x,y) - change as desired
% Example 1: combination of sines
u0 = @(x,y) sin(pi*x/Lx).*sin(2*pi*y/Ly) + 0.5*sin(2*pi*x/Lx).*sin(pi*y/Ly);

% If you'd like a different u0, replace the line above with your function handle.
% For instance: u0 = @(x,y) exp(-50*((x-0.5*Lx).^2 + (y-0.5*Ly).^2));  % gaussian

% Precompute sine basis on spatial grid (for reconstruction)
sinx = zeros(M, Nx);
siny = zeros(N, Ny);
for m = 1:M
    sinx(m,:) = sin(m*pi*x/Lx);
end
for n = 1:N
    siny(n,:) = sin(n*pi*y/Ly);
end

% Compute coefficients A(m,n) numerically using double trapezoidal integration
% A_mn = (4/(Lx*Ly)) * \int_0^{Lx}\int_0^{Ly} u0(x,y) sin(m*pi*x/Lx) sin(n*pi*y/Ly) dy dx
fprintf('Computing coefficients A_{mn} (M=%d, N=%d)...\n', M, N);
% Create fine quadrature grid for integration (can be coarser if u0 smooth)
qx = linspace(0,Lx,2*Nx); qy = linspace(0,Ly,2*Ny);
[QX,QY] = meshgrid(qx,qy);
U0vals = u0(QX,QY);

A = zeros(M,N);
wx = (Lx/(length(qx)-1)); wy = (Ly/(length(qy)-1)); % trapezoid spacing factors (we'll use trapz)
% We'll use trapz for each mode
for m = 1:M
    sx = sin(m*pi*QX/Lx);  % same size as QX
    for n = 1:N
        sy = sin(n*pi*QY/Ly);
        integrand = U0vals .* sx .* sy;
        % double integral with trapz: first integrate w.r.t y (columns), then x (rows)
        int_y = trapz(qy, integrand, 1);   % integrate along rows -> column vector over x
        I = trapz(qx, int_y);              % final scalar
        A(m,n) = (4/(Lx*Ly)) * I;
    end
end

% Reconstruct solution at requested times
for tidx = 1:length(t_plot)
    t = t_plot(tidx);
    U = zeros(Ny, Nx);
    % sum over modes
    for m = 1:M
        for n = 1:N
            lambda = k * pi^2 * (m^2/Lx^2 + n^2/Ly^2);
            coeff = A(m,n) * exp(-lambda * t);
            % outer product sin(m*pi*x/Lx)*sin(n*pi*y/Ly)
            % sinx(m,:) is 1 x Nx, siny(n,:)' is Ny x 1
            U = U + coeff * (siny(n,:)' * sinx(m,:));
        end
    end

    % Plot
    figure(tidx);
    surf(X,Y,U, 'EdgeColor','none');
    xlabel('x'); ylabel('y'); zlabel('u(x,y,t)');
    title(sprintf('2D Heat (series) solution at t = %.4f', t));
    view(35,30);
    axis([0 Lx 0 Ly min(U(:)) max(U(:))]);
    colorbar;
    drawnow;
end

% Optional: animation
figure;
for tidx = 1:100
    t = (tidx-1) * 0.02;  % change final time/step as desired
    U = zeros(Ny, Nx);
    for m = 1:M
        for n = 1:N
            lambda = k * pi^2 * (m^2/Lx^2 + n^2/Ly^2);
            coeff = A(m,n) * exp(-lambda * t);
            U = U + coeff * (siny(n,:)' * sinx(m,:));
        end
    end
    surf(X,Y,U,'EdgeColor','none'); title(sprintf('t = %.3f',t));
    xlabel('x'); ylabel('y'); zlabel('u'); view(35,30); colorbar;
    axis([0 Lx 0 Ly -1 1]); drawnow;
end
