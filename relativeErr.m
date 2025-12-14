clc; clear; close all;

% Time values
t = [0.05 0.10 0.15 0.25];

% Relative error data
max_err = [0.00110 0.00195 0.00245 0.00382];
min_err = [0.00002 0.00003 0.00005 0.00007];
avg_err = [0.00048 0.00079 0.00112 0.00194];

% Plot
figure; hold on;

plot(t, max_err, 'r-o', 'LineWidth', 2, 'MarkerSize', 8);
plot(t, min_err, 'b-s', 'LineWidth', 2, 'MarkerSize', 8);
plot(t, avg_err, 'g-^', 'LineWidth', 2, 'MarkerSize', 8);

title('Relative Error Analysis of FTCS Method', 'FontSize', 16);
xlabel('Time (t)', 'FontSize', 14);
ylabel('Relative Error', 'FontSize', 14);

legend('Maximum Relative Error', 'Minimum Relative Error', 'Average Relative Error', ...
       'Location', 'northwest');

grid on;
hold off;
