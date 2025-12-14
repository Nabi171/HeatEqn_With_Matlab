% Parameters
k = 1;

% Spatial and time grids
x = linspace(0, 1, 200);
t = linspace(0, 0.25, 200);
[X, T] = meshgrid(x, t);

% Analytical solution
U = sin(pi * X) .* exp(-k * (pi^2) * T);

% Create 3D surface plot
figure;
s = surf(X, T, U);
shading interp;
colormap jet;
colorbar;
hold on;

% ==== ADD LEGEND TRICK ====
% Create invisible plot handles for legend entries
p1 = plot3(nan, nan, nan, 'Color', 'b', 'LineWidth', 2);
p2 = plot3(nan, nan, nan, 'Color', 'r', 'LineWidth', 2);
p3 = plot3(nan, nan, nan, 'Color', 'g', 'LineWidth', 2);
p4 = plot3(nan, nan, nan, 'Color', 'm', 'LineWidth', 2);

legend([p1 p2 p3 p4], ...
    'High temperature region', ...
    'Medium-high region', ...
    'Medium-low region', ...
    'Low temperature region', ...
    'Location', 'northeastoutside');

% Labels & title
xlabel('Position (x)', 'FontSize', 14);
ylabel('Time (t)', 'FontSize', 14);
zlabel('Temperature u(x,t)', 'FontSize', 14);
title('3D Surface Plot of Analytical Heat Equation Solution', 'FontSize', 16);

grid on;
view(45, 30);
hold off;
