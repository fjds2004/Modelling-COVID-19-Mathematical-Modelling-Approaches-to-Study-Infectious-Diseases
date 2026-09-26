clear; close all; clc;

% Population
N = 1000;           


beta = 0.5;         % Transmission rate (per day).
sigma = 1/5;      % Rate of progression from E to I.
gamma = 0.2;      % Recovery rate.

% Calculate R0
R0 = beta/gamma;    


% Average time in exposed state
latent_period = 1/sigma;
% Average time in infected state
infectious_period = 1/gamma;   



I0 = 1;             % Initial infected
E0 = 0;             % Initial exposed (none at start)
R0_init = 0;        % Initial recovered
S0 = N - E0 - I0 - R0_init;  % Initial susceptible

% Time parameters
t_span = [0 200];   
t_eval = linspace(0, 200, 1000);


fprintf('SEIR Model Parameters:\n');
fprintf('===================================\n');
fprintf('Population (N): %d\n', N);
fprintf('Transmission rate (β): %.3f per day\n', beta);
fprintf('Latent period: %.1f days\n', latent_period);
fprintf('Infectious period: %.1f days\n', infectious_period);
fprintf('Progression rate (σ): %.3f per day\n', sigma);
fprintf('Recovery rate (γ): %.2f per day\n', gamma);
fprintf('Basic reproduction number (R₀): %.2f\n', R0);
fprintf('Initial conditions: S=%d, E=%d, I=%d, R=%d\n', ...
    S0, E0, I0, R0_init);
fprintf('===================================\n\n');

%% Solve the ODE System
y0 = [S0; E0; I0; R0_init];

% Use ODE45 solver
[t, y] = ode45(@(t,y) ode_SEIR(t, y, beta, sigma, gamma, ...
    N), t_eval, y0);

% Compartments
S = y(:,1);  % Susceptible
E = y(:,2);  % Exposed
I = y(:,3);  % Infected
R = y(:,4);  % Recovered

%% Find key metrics
[peak_exposed, peak_E_idx] = max(E);
peak_E_time = t(peak_E_idx);

[peak_infected, peak_I_idx] = max(I);
peak_I_time = t(peak_I_idx);

final_recovered = R(end);
attack_rate = (final_recovered / N) * 100;

fprintf('Epidemic Outcomes:\n');
fprintf('===================================\n');
fprintf(['Peak exposed: %.0f individuals (%.1f%% of ' ...
    'population)\n'],peak_exposed, (peak_exposed/N)*100);
fprintf('Time to peak exposed: %.1f days\n', peak_E_time);
fprintf(['Peak infected: %.0f individuals (%.1f%% of ' ...
    '    population)\n'],peak_infected, (peak_infected/N) ...
    *100);
fprintf('Time to peak infected: %.1f days\n', peak_I_time);
fprintf('Time delay (E peak to I peak): %.1f days\n', ...
    peak_I_time - peak_E_time);
fprintf('Final attack rate: %.1f%%\n', attack_rate);
fprintf('Duration of epidemic: ~%.0f days\n', ...
    t(find(I > 1, 1, 'last')));
fprintf('===================================\n');

%% Create figure
figure('Position', [100 100 900 600]);

% Main plot: all four compartments
hold on;
plot(t, S, 'b-', 'LineWidth', 2.5, 'DisplayName', ...
    'Susceptible (S)');
plot(t, E, 'm-', 'LineWidth', 2.5, 'DisplayName', ...
    'Exposed (E)');
plot(t, I, 'r-', 'LineWidth', 2.5, 'DisplayName', ...
    'Infected (I)');
plot(t, R, 'g-', 'LineWidth', 2.5, 'DisplayName', ...
    'Recovered (R)');

% Mark peak exposed and infected
plot(peak_E_time, peak_exposed, 'black*', ...
    'MarkerSize', 15, 'LineWidth', 2, ...
     'HandleVisibility', 'off');
plot(peak_I_time, peak_infected, 'black*', ...
    'MarkerSize', 15, 'LineWidth', 2, ...
     'HandleVisibility', 'off');

% Formatting
xlabel('Time (days)', 'FontSize', 14, ...
    'FontWeight', 'bold');
ylabel('Number of Individuals', 'FontSize', 14, ...
    'FontWeight', 'bold');
title(sprintf(['SEIR Epidemic Model: R_0 = %.2f, ' ...
    '    Latent Period = %.1f days'], ...
      R0, latent_period), 'FontSize', 16, ...
      'FontWeight', 'bold');
legend('Location', 'east', 'FontSize', 12);
grid on;
set(gca, 'FontSize', 12);
xlim([0 200]);
ylim([0 N]);

% Add text box with key information
annotation('textbox', [0.645, 0.63, 0.25, 0.15], ...
    'String', {sprintf('R_0 = %.2f', R0), ...
               sprintf('Latent period = %.1f days', latent_period), ...
               sprintf('Peak (E): %.0f at T=%.1f days', peak_exposed, peak_E_time), ...
               sprintf('Peak (I): %.0f at T=%.1f days', peak_infected, peak_I_time), ...
               sprintf('Attack rate: %.1f%%', attack_rate)}, ...
    'FontSize', 10, ...
    'BackgroundColor', 'white', ...
    'EdgeColor', 'black', ...
    'LineWidth', 0.5);

hold off;

