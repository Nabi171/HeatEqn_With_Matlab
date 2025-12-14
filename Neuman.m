% 1D Heat equation - Neumann boundary condition
% u_t = k u_xx,  u_x(0,t)=u_x(L,t)=0

clear; close all; clc;

% Parameters
L = 1;        % domain length
k = 1;        % thermal diffusivity
N = 40;       % number of cosine modes

% Spatial grid
x = linspace(0,L,400);

% Time values
t_vals = [0, 0.01, 0.05, 0.2];

% Initial condition (change freely)
u0 = @(x) cos(pi*x/L) + 0.5*cos(3*pi*x/L) + 1;

% Compute coefficients
A0 = (1/L) * trapz(x, u0(x));
A  = zeros(1,N);

for n = 1:N
    integrand = u0(x).*cos(n*pi*x/L);
    A(n) = (2/L)*trapz(x, integrand);
end

% Plot solution at different times
figure;
hold on; grid on;

for t = t_vals
    u = A0*ones(size(x));
    for n = 1:N
        u = u + A(n)*cos(n*pi*x/L).*exp(-k*(n*pi/L)^2*t);
    end
    plot(x, u, 'LineWidth', 2, 'DisplayName', sprintf('t = %.2f',t));
end

xlabel('x');
ylabel('u(x,t)');
title('1D Heat Equation (Neumann Boundary Condition)');
legend;
