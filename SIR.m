clear; close all; clc;

% Total population size
N=1000;           

% Transmission rate (per day)
beta=0.35;

% Recovery rate (per day)
gamma=0.14;

% Basic reproduction number
R0=beta/gamma;   


% starting conditions
I0 = 1;             
R0_start = 0;        
S0 = N - I0 - R0_start;   

% Time parameters
t_span = [0 200];   % 200 days

t_eval = linspace(0, 200, 1000);

%% parameters display
fprintf('SIR Model Parameters:\n');
fprintf('----------------------------\n');
fprintf('Population (N): %d\n', N);
fprintf('Transmission rate (β): %.2f per day\n', beta);
fprintf('Recovery rate (γ): %.2f per day\n', gamma);
fprintf('Infectious period: %.1f days\n', 1/gamma);
fprintf('Basic reproduction number (R0): %.2f\n', R0);
fprintf('Initial infected: %d\n', I0);
fprintf('----------------------------\n\n');

%% Solve the ODE System
y0 = [S0; I0; R0_start];

% Use ODE45 solver
[t, y] = ode45(@(t,y) ode_SIR(t,y,beta,gamma,N), t_eval, y0);

% Extract compartments
S = y(:,1);  % Susceptible
I = y(:,2);  % Infected
R = y(:,3);  % Recovered


%% Find key metrics
[peak_infected, peak_idx] = max(I);
peak_time = t(peak_idx);
final_recovered = R(end);
attack_rate = (final_recovered / N) * 100;

fprintf('Epidemic Outcomes:\n');
fprintf('----------------------------\n');
fprintf('Peak infected: %.0f individuals (%.1f%% of population)\n', ...
        peak_infected, (peak_infected/N)*100);
fprintf('Time to peak: %.1f days\n', peak_time);
fprintf('Final attack rate: %.1f%%\n', attack_rate);
fprintf('Duration of epidemic: ~%.0f days\n', t(find(I > 1, 1, 'last')));
fprintf('----------------------------\n');

%% Create figure
figure('Position', [100 100 900 600]);

% Main plot: all three compartments
hold on;
plot(t, S, 'blue', 'LineWidth', 2.5, 'DisplayName', 'Susceptible (S)');
plot(t, I, 'red', 'LineWidth', 2.5, 'DisplayName', 'Infected (I)');
plot(t, R, 'green', 'LineWidth', 2.5, 'DisplayName', 'Recovered (R)');

% Mark peak infection
plot(peak_time, peak_infected, 'Black*', 'MarkerSize', 15, ...
    'LineWidth', 2,'DisplayName', sprintf('Peak: %.0f at T= %.1f days', ...
    peak_infected, peak_time));

% Formatting
xlabel('Time (days)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Number of Individuals', 'FontSize', 14, 'FontWeight', 'bold');
title(sprintf('SIR Epidemic Model: R_0 = %.2f', R0), 'FontSize', ...
    16, 'FontWeight', 'bold');
legend('Location', 'east', 'FontSize', 12);
grid on;
set(gca, 'FontSize', 12);
xlim([0 200]);
ylim([0 N]);

hold off;

%% Save figure
print('fig_SIR_epidemic_curve', '-dpng', '-r300');
fprintf('\nFigure saved as: fig_SIR_epidemic_curve.png\n');

title('SIR Model of COVID-19 Epidemic Dynamics (R_0 = 2.5)', ...
    'FontSize', 16, 'FontWeight', 'bold');