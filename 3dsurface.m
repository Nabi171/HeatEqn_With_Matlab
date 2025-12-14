% Parameters
k = 1;

% Spatial and time grids
x = linspace(0, 1, 200);
t = linspace(0, 3, 200);

[X, T] = meshgrid(x, t);

% Analytical solution u(x,t)
U = sin(pi * X) .* exp(-k * (pi^2) * T);

% Plot 3D surface
figure;
surf(X, T, U);

% Labels and title
xlabel('Position x', 'FontSize', 12);
ylabel('Time t', 'FontSize', 12);
zlabel('Temperature u(x,t)', 'FontSize', 12);

title('3D Surface Plot of Analytical Heat Equation Solution', 'FontSize', 14);

% Improve visual quality
shading interp;         % Smooth the color transitions
colormap jet;           % Color pattern
colorbar;               % Add color scale
view(45, 30);           % Better viewing angle
grid on;
