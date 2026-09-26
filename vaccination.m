clear; close all; clc;

%% Model Parameters
N = 1000;
gamma = 1/10;       % 10 day infectious period
delta = 1/5;        % 5 day incubation period
R0 = 2.5;
beta = R0 * gamma;  % Transmission rate


T = 365;            % Simulation (days/year)
t_eval = linspace(0, T, 1000);

%% SEIR Model (no vaccination)
function dydt = SEIR_equations(t, y, beta, delta, gamma, N)
    S = y(1);
    E = y(2);
    I = y(3);
    R = y(4);
    
    dS = -beta * S * I / N;
    dE = beta * S * I / N - delta * E;
    dI = delta * E - gamma * I;
    dR = gamma * I;
    
    dydt = [dS; dE; dI; dR];
end

%% Pre-Vaccination (reduce S(0) before epidemic starts)

vaccination_levels = [0, 0.20, 0.40, 0.60];
colors = {[1 0 0], [1 0.5 0], [0 0.6 0], [0 0 1]};
linestyles = {'-', '-', '-', '-'};

figure('Position', [100 100 800 500]);
hold on;
for k = 1:length(vaccination_levels)
    v = vaccination_levels(k);
    
    
    % Initial conditions: move fraction v from S to R
    S0 = N * (1 - v) - 1;
    E0 = 0;
    I0 = 1;
    R0_init = N * v;
    
    y0 = [S0; E0; I0; R0_init];
    
    [t, y] = ode45(@(t,y) SEIR_equations(t, y, beta, delta, gamma, N), t_eval, y0);
    
    I_curve = y(:,3);
    [peak_I, peak_idx] = max(I_curve);
    peak_t = t(peak_idx);
    
    Rv = (1 - v) * R0;
    label = sprintf('v = %.2f (R_v = %.2f)', v, Rv);
    
    plot(t, I_curve, 'Color', colors{k}, 'LineStyle', linestyles{k}, ...
         'LineWidth', 3, 'DisplayName', label);
    
    % Mark peak with asterisk
    if peak_I > 5
        plot(peak_t, peak_I, '*', 'Color', colors{k}, 'MarkerSize', 15, ...
             'LineWidth', 2, 'HandleVisibility', 'off');
    end

end

% Formatt
xlabel('Time (days)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Infected Individuals', 'FontSize', 14, 'FontWeight', 'bold');
title(sprintf('Effect of Vaccination on SEIR Epidemic Dynamics (R_0 = 2.50)', R0), ...
      'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'northeast', 'FontSize', 11);
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
xlim([0 T]);
ylim([0 165]);
hold off;


