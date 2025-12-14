% Parameters
k = 1;
x = linspace(0, 1, 500);

% Time levels you want to show
times = [0, 1, 2, 3];

% Colors for each individual figure
colors = {'b', [1 0.5 0], 'g', 'r'};   % blue, orange, green, red

for i = 1:length(times)
    
    t = times(i);
    u = sin(pi * x) .* exp(-k * (pi^2) * t);
    
    figure;
    plot(x, u, 'Color', colors{i}, 'LineWidth', 2);
    
    title(['Heat Equation Analytical Solution at t = ' num2str(t)], 'FontSize', 14);
    xlabel('Position (x)', 'FontSize', 12);
    ylabel('Temperature u(x,t)', 'FontSize', 12);
    
    grid on;
    
end
