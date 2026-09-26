clear; close all; clc;

%% Model Parameters
N = 1000;           
gamma = 0.2;        
beta = 0.6;         
R0 = beta/gamma;

% Calculate endemic equilibrium
I_star = N * (1 - 1/R0);

% Initial conditions
I0 = 100;           
S0 = N - I0;        

fprintf('SIS Model Parameters (R0 > 1):\n');
fprintf('================================\n');
fprintf('Beta (β): %.2f per day\n', beta);
fprintf('Gamma (γ): %.2f per day\n', gamma);
fprintf('R0: %.2f\n', R0);
fprintf('Initial infected: %d (%.1f%%)\n', I0, 100*I0/N);
fprintf('Endemic equilibrium I*: %.0f (%.1f%% of population)\n', ...
    I_star, 100*I_star/N);
fprintf('================================\n\n');

%% Solve the ODE
t_span = [0 200];
t_eval = linspace(0, 200, 1000);

[t, y] = ode45(@(t,y) SIS_equations(t, y, beta, gamma, N), t_eval, [S0; I0]);

S = y(:,1);  
I = y(:,2);  

%% Plot

figure('Position', [100 100 800 600]);

% Plot infected population
plot(t, I, 'b-', 'LineWidth', 3);
hold on;

% Horizontal line endemic equilibrium
plot([0 200], [I_star I_star], 'r--', 'LineWidth', 2, ...
    'DisplayName', sprintf('Endemic Equilibrium'));

% Format
xlabel('Time (days)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Infected Population', 'FontSize', 14, 'FontWeight', 'bold');
title(sprintf('SIS Model: Endemic Equilibrium (R_0 = %.2f > 1)', R0), ...
    'FontSize', 16, 'FontWeight', 'bold');
legend('Infected I(t)', sprintf('Endemic Equilibrium', I_star), ...
    'Location', 'southeast', 'FontSize', 11);
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
xlim([0 200]);
ylim([0 N*0.75]);

% Add text 
text(100, I_star*1.3, {sprintf('R_0 = %.2f > 1', R0), ...
                        'Disease persists', sprintf('I(t) → I^* = %.0f', ...
                        I_star)},FontSize', 12, 'FontWeight', ...
                        'bold', 'BackgroundColor', 'white','EdgeColor', ...
                        'black', 'LineWidth', 1.5,'HorizontalAlignment', ...
                        'center');

hold off;

%% Save figure
print('fig_SIS_R0_greater_than_1', '-dpng', '-r300');
fprintf('Figure saved as: fig_SIS_R0_greater_than_1.png\n');

fprintf(['Infection reaches endemic equilibrium at I* = %.0f ' ...
    '        (%.1f%% of population)\n'], I_star, 100*I_star/N);
fprintf('\nFigure saved as: fig_SIS_epidemic_R>1.png\n');


function dydt = SIS_equations(t, y, beta, gamma, N)
    % SIS model differential equations
    
    S = y(1);  % Susceptible
    I = y(2);  % Infected
    
    % Force of infection
    lambda = beta * I / N;
    
    % Differential equations
    dS = -lambda * S + gamma * I;  % Susceptibles decrease by infection
    dI = lambda * S - gamma * I;    % Infected increase by infection
    
    dydt = [dS; dI];
end