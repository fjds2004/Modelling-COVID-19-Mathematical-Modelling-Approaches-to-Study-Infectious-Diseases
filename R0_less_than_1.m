clear; close all; clc;

%% Model Parameters
N = 1000;           
gamma = 0.2;        
beta = 0.15;        
R0 = beta/gamma;    

% Initial conditions
I0 = 100;           
S0 = N - I0;        

fprintf('SIS Model Parameters (R0 < 1):\n');
fprintf('================================\n');
fprintf('Beta (β): %.2f per day\n', beta);
fprintf('Gamma (γ): %.2f per day\n', gamma);
fprintf('R0: %.2f\n', R0);
fprintf('Initial infected: %d (%.1f%%)\n', I0, 100*I0/N);
fprintf('================================\n\n');

%% Solve the ODE
t_span = [0 200];
t_eval = linspace(0, 200, 1000);

[t, y] = ode45(@(t,y) SIS_equations(t, y, beta, gamma, N), t_eval, [S0; I0]);

S = y(:,1);  
I = y(:,2);  

%% Create Figure
figure('Position', [100 100 800 600]);

% Plot infected population
plot(t, I, 'b-', 'LineWidth', 3);
hold on;

% Add horizontal line at I = 0
plot([0 200], [0 0], 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% Format
xlabel('Time (days)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Infected Population', 'FontSize', 14, 'FontWeight', 'bold');
title(sprintf('SIS Model: Disease Dies Out (R_0 = %.2f < 1)', R0), ...
      'FontSize', 16, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
xlim([0 200]);
ylim([0 I0*1.1]);

% text annotation
text(100, I0*0.7, {sprintf('R_0 = %.2f < 1', R0), ...
                   'Disease cannot persist', ...
                   'I(t) → 0 as t → ∞'}, ...
     'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', 'white', ...
     'EdgeColor', 'black', ...
     'LineWidth', 1.5, ...
     'HorizontalAlignment', 'center');

hold off;

%% Save figure
print('fig_SIS_R0_less_than_1', '-dpng', '-r300');
fprintf('Figure saved as: fig_SIS_R0_less_than_1.png\n');
fprintf('Infection dies out to zero (disease-free equilibrium is stable)\n');
fprintf('\nFigure saved as: fig_SIS_epidemic_R<1.png\n');


function dydt = SIS_equations(t, y, beta, gamma, N)
    
    S = y(1);  
    I = y(2);  
    
    % Force of infection
    lambda = beta * I / N;
    
    % Differential equations
    dS = -lambda * S + gamma * I;  
    dI = lambda * S - gamma * I;    
    
    dydt = [dS; dI];
end
