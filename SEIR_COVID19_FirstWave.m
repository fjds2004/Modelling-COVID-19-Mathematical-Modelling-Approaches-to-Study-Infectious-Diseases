function SEIR_COVID19_FirstWave()
close all; clc;

%% Population
N = 67100000;       % UK population 
%% Deriving Parameters from Reported Data


% Incubation period: estimated at 5.2 days
sigma = 1 / 5.2;

% (2) Infectious period
gamma = 1 / 10;

% (3) Basic reproduction number early UK estimates placed R0
%     between 2.5 and 3.5
R0 = 3.0;
beta = R0 * gamma;


I0 = 20000;                       % Estimated true infections
E0 = round(I0 * gamma / sigma);   % Derived exposed count
R0_init = 0;
S0 = N - E0 - I0 - R0_init;


% Calculate derived quantities
latent_period     = 1 / sigma;
infectious_period = 1 / gamma;


%% Time parameters (first wave = 219 days)
T_days = 219;
t_eval = linspace(0, T_days, 2000);

fprintf('SEIR Model for COVID-19 (First Wave):\n');
fprintf('===================================\n');
fprintf('Population (N): %d\n', N);
fprintf('Transmission rate (beta): %.4f per day\n', beta);
fprintf('Latent period (1/sigma): %.1f days\n', latent_period);
fprintf('Infectious period (1/gamma): %.1f days\n', infectious_period);
fprintf('Progression rate (sigma): %.4f per day\n', sigma);
fprintf('Recovery rate (gamma): %.2f per day\n', gamma);
fprintf('Basic reproduction number (R0): %.2f\n', R0);
fprintf('Derived initial conditions:\n');
fprintf('  S(0) = %d\n', S0);
fprintf('  E(0) = %d  (from E/I = gamma/sigma)\n', E0);
fprintf('  I(0) = %d  (estimated true infections)\n', I0);
fprintf('  R(0) = %d\n', R0_init);
fprintf('===================================\n\n');


%% Solve the ODE System
y0 = [S0; E0; I0; R0_init];

[t, y] = ode45(@(t,y) ode_SEIR(t, y, beta, sigma, gamma, ...
    N), t_eval, y0);

% Compartments
S = y(:,1);   % Susceptible
E = y(:,2);   % Exposed
I = y(:,3);   % Infected
R = y(:,4);   % Recovered

%% Find key metrics
[peak_exposed, peak_E_idx] = max(E);
peak_E_time = t(peak_E_idx);

[peak_infected, peak_I_idx] = max(I);
peak_I_time = t(peak_I_idx);

final_recovered = R(end);
attack_rate = (final_recovered / N) * 100;

fprintf('Epidemic Outcomes (First Wave):\n');
fprintf('===================================\n');
fprintf('Peak exposed: %.0f individuals\n', peak_exposed);
fprintf('Time to peak exposed: %.1f days\n', peak_E_time);
fprintf('Peak infected: %.0f individuals\n', peak_infected);
fprintf('Time to peak infected: %.1f days\n', peak_I_time);
fprintf('Time delay (E peak to I peak): %.1f days\n', ...
    peak_I_time - peak_E_time);
fprintf('Final attack rate: %.1f%%\n', attack_rate);
fprintf('R0 = beta/gamma = %.2f\n', R0);
fprintf('===================================\n');

%% Create figure
figure('Position', [100 100 900 600]);

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
title(sprintf(['SEIR Model: COVID-19 First Wave ' ...
    '(R_0 = %.2f, Latent Period = %.1f days)'], ...
    R0, latent_period), 'FontSize', 16, ...
    'FontWeight', 'bold');
legend('Location', 'east', 'FontSize', 12);
grid on;
set(gca, 'FontSize', 12);
xlim([0 T_days]);
ylim([0 N]);

% Add text box with key information
annotation('textbox', [0.594, 0.55, 0.3, 0.22], ...
    'String', {sprintf('R_0 = %.2f', R0), ...
               sprintf('Latent period = %.1f days', ...
               latent_period), ...
               sprintf('Infectious period = %.1f days', ...
               infectious_period), ...
               sprintf('E(0) = %d,  I(0) = %d', E0, I0), ...
               sprintf('Peak (E): %.0f at T=%.1f days', ...
               peak_exposed, peak_E_time), ...
               sprintf('Peak (I): %.0f at T=%.1f days', ...
               peak_infected, peak_I_time), ...
               sprintf('Attack rate: %.1f%%', attack_rate)}, ...
    'FontSize', 10, ...
    'BackgroundColor', 'white', ...
    'EdgeColor', 'black', ...
    'LineWidth', 0.5);
hold off;

%% Save figure
print('fig_SEIR_COVID19_first_wave', '-dpng', '-r300');
fprintf('\nFigure saved as: fig_SEIR_COVID19_first_wave.png\n');

end

%% ODE function
function dydt = ode_SEIR(~, y, beta, sigma, gamma, N)
    S = y(1);
    E = y(2);
    I = y(3);
    R = y(4);

    dS = -beta * S * I / N;
    dE =  beta * S * I / N - sigma * E;
    dI =  sigma * E - gamma * I;
    dR =  gamma * I;

    dydt = [dS; dE; dI; dR];
end