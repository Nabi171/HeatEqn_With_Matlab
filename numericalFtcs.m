% Numerical FTCS solution of the 1D Heat Equation
k = 1;               % diffusion coefficient
L = 1;               % rod length
Nx = 100;            % number of spatial intervals
dx = L / Nx;         % spatial step
x = linspace(0, L, Nx+1);

dt = 0.0001;         % time step (satisfies stability r <= 0.5)
r = k * dt / dx^2;   % Fourier number

fprintf("Fourier number r = %f\n", r);

Nt = 1500;           % total time steps
t_values = [0, 0.05, 0.10, 0.15];  % times to plot
t_steps = round(t_values / dt);    % convert time to steps

% Initial condition
u = sin(pi * x);

% Storage for plotting states
U_plot = zeros(length(t_values), length(x));
U_plot(1, :) = u;

% Time loop (FTCS update)
for n = 1:Nt
    u_new = u;
    
    for i = 2:Nx    % interior points
        u_new(i) = u(i) + r * (u(i+1) - 2*u(i) + u(i-1));
    end
    
    % Apply boundary conditions
    u_new(1) = 0;
    u_new(end) = 0;
    
    u = u_new;

    % Save solutions at specified times
    for ktime = 2:length(t_values)
        if n == t_steps(ktime)
            U_plot(ktime, :) = u;
        end
    end
end

% === Plot FTCS results ===
figure;
hold on;

plot(x, U_plot(1, :), 'b', 'LineWidth', 2);
plot(x, U_plot(2, :), 'r', 'LineWidth', 2);
plot(x, U_plot(3, :), 'g', 'LineWidth', 2);
plot(x, U_plot(4, :), 'm', 'LineWidth', 2);

xlabel('x', 'FontSize', 14);
ylabel('u(x,t)', 'FontSize', 14);

title('Numerical FTCS Solution of the Heat Equation', 'FontSize', 16);

legend('t = 0', 't = 0.05', 't = 0.10', 't = 0.15', ...
       'Location', 'northeastoutside');

grid on;
hold off;
