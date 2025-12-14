k = 1;
x = linspace(0, 1, 500);

% Time values for separate colors
times = [0, 0.05, 0.10, 0.15];
colors = {'b', 'r', 'g', 'm'}; % blue, red, green, magenta

figure;
hold on;

for i = 1:length(times)
    t = times(i);
    u = sin(pi .* x) .* exp(-k * (pi^2) * t);
    plot(x, u, 'Color', colors{i}, 'LineWidth', 2, 'DisplayName', ['t = ' num2str(t)]);
end

xlabel('x', 'FontSize', 12);
ylabel('u(x,t)', 'FontSize', 12);
title('Heat Equation Solution for Different Time Values', 'FontSize', 14);
grid on;
legend show;
hold off;
