% Parameters
k = 1;

% Spatial and time grids
x = linspace(0, 1, 200);
t = linspace(0, 0.25, 200); % adjust as needed
[X, T] = meshgrid(x, t);

% Analytical solution u(x,t)
U = sin(pi * X) .* exp(-k * (pi^2) * T);

% Create 3D surface plot
figure;
surf(X, T, U);

% Labels and title
xlabel('Position (x)', 'FontSize', 14);
ylabel('Time (t)', 'FontSize', 14);
zlabel('Temperature u(x,t)', 'FontSize', 14);
title('3D Surface Plot of the Analytical Heat Equation Solution', 'FontSize', 16);

% Improve appearance
shading interp;       % smooth color transitions
colormap jet;         % colorful heatmap colormap
colorbar;             % add colorbar on side
grid on;              % show grid lines
view(45, 30);         % set 3D viewing angle
