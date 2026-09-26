%% Extended SEIR Model Simulation (with Vaccination and Deceased classes)
%  Building on He et al. (2020) framework with V and D compartments
%  Frankie Docking-Smith - Dissertation
%  ---------------------------------------------------------------

clear; clc; close all;

%% ===================== PARAMETERS =====================

% --- Population ---
N0 = 1000;          % Initial total population

% --- Transmission rates ---
beta1 = 0.3;        % Contact/infection rate from I1 (undetected)
beta2 = 0.1;        % Contact/infection rate from I2 (detected, lower due to intervention)
nu    = 0.05;       % Transmission probability from exposed individuals

% --- Progression rates ---
theta1 = 0.1;       % Rate from E to I1 (undetected infectious)
theta2 = 0.05;      % Rate from E to I2 (detected infectious)

% --- Recovery rates ---
gamma1 = 0.07;      % Recovery rate for I1
gamma2 = 0.05;      % Recovery rate for I2
phi    = 0.1;       % Recovery rate for hospitalised (H)

% --- Hospitalisation ---
varphi = 0.05;      % Rate of hospitalisation from I2

% --- Quarantine ---
q1 = 0.02;          % Rate from Q back to S
q2 = 0.01;          % Rate from S to Q
lambda_q = 0.03;    % Rate from Q to I2

% --- Waning immunity ---
alpha = 0.005;      % Rate at which recovered lose immunity (R -> S)

% --- External input ---
K = 2;              % Imported cases per day

% --- Vaccination (NEW) ---
rho = 0.005;        % Vaccination rate (S -> V)

% --- Death rates (NEW) ---
mu1 = 0.005;        % Death rate for I1 (no intervention)
mu2 = 0.003;        % Death rate for I2 (with intervention)
mu3 = 0.002;        % Death rate for H (hospitalised)

%% ===================== INITIAL CONDITIONS =====================

S0  = N0 - 20;     % Susceptible
E0  = 10;           % Exposed
I1_0 = 5;           % Infectious without intervention
I2_0 = 5;           % Infectious with intervention
R0_init = 0;        % Recovered
Q0  = 0;            % Quarantined
H0  = 0;            % Hospitalised
V0  = 0;            % Vaccinated
D0  = 0;            % Deceased

y0 = [S0, E0, I1_0, I2_0, R0_init, Q0, H0, V0, D0];

%% ===================== SOLVE ODEs =====================

tspan = [0 365];    % Simulate for 1 year

% Pack parameters into a struct for clarity
params.beta1 = beta1;   params.beta2 = beta2;   params.nu = nu;
params.theta1 = theta1; params.theta2 = theta2;
params.gamma1 = gamma1; params.gamma2 = gamma2;  params.phi = phi;
params.varphi = varphi;
params.q1 = q1;         params.q2 = q2;         params.lambda_q = lambda_q;
params.alpha = alpha;    params.K = K;
params.rho = rho;
params.mu1 = mu1;       params.mu2 = mu2;       params.mu3 = mu3;

[t, y] = ode45(@(t, y) extended_seir(t, y, params), tspan, y0);

% Extract variables
S  = y(:,1);  E  = y(:,2);  I1 = y(:,3);  I2 = y(:,4);
R  = y(:,5);  Q  = y(:,6);  H  = y(:,7);  V  = y(:,8);  D = y(:,9);

% Total living population
N_living = S + E + I1 + I2 + R + Q + H + V;

% Total infected (I1 + I2)
I_total = I1 + I2;

%% ===================== FIND KEY METRICS =====================

% Peak total infection
[peak_I, idx_peak] = max(I_total);
t_peak = t(idx_peak);

% Final values
fprintf('=============== SIMULATION RESULTS ===============\n');
fprintf('Peak total infected: %.0f at T = %.1f days\n', peak_I, t_peak);
fprintf('Final recovered:     %.0f\n', R(end));
fprintf('Final vaccinated:    %.0f\n', V(end));
fprintf('Final deceased:      %.0f\n', D(end));
fprintf('Final living pop:    %.0f (started at %d)\n', N_living(end), N0);
fprintf('==================================================\n');

%% ===================== PLOTTING =====================

% --- Figure 1: All compartments ---
figure('Position', [100, 100, 900, 600]);

subplot(2,2,1)
plot(t, S, 'b-', 'LineWidth', 1.5); hold on;
plot(t, E, 'Color', [1 0.5 0], 'LineWidth', 1.5);
plot(t, R, 'g-', 'LineWidth', 1.5);
xlabel('Time (days)'); ylabel('Number of Individuals');
title('Susceptible, Exposed, and Recovered');
legend('S', 'E', 'R', 'Location', 'best');
grid on;

subplot(2,2,2)
plot(t, I1, 'r-', 'LineWidth', 1.5); hold on;
plot(t, I2, 'm-', 'LineWidth', 1.5);
plot(t, I_total, 'k--', 'LineWidth', 1.5);
plot(t_peak, peak_I, 'k*', 'MarkerSize', 12, 'LineWidth', 2);
xlabel('Time (days)'); ylabel('Number of Individuals');
title('Infectious Populations');
legend('I_1 (no intervention)', 'I_2 (with intervention)', ...
       'Total infected', sprintf('Peak: %.0f at T=%.1f', peak_I, t_peak), ...
       'Location', 'best');
grid on;

subplot(2,2,3)
plot(t, Q, 'Color', [0.5 0 0.5], 'LineWidth', 1.5); hold on;
plot(t, H, 'Color', [0 0.5 0.5], 'LineWidth', 1.5);
plot(t, V, 'Color', [0 0.7 0], 'LineWidth', 1.5);
xlabel('Time (days)'); ylabel('Number of Individuals');
title('Quarantined, Hospitalised, and Vaccinated');
legend('Q', 'H', 'V', 'Location', 'best');
grid on;

subplot(2,2,4)
plot(t, D, 'k-', 'LineWidth', 2); hold on;
plot(t, N_living, 'b--', 'LineWidth', 1.5);
xlabel('Time (days)'); ylabel('Number of Individuals');
title('Deceased and Living Population');
legend('Deceased (D)', 'Living population N(t)', 'Location', 'best');
grid on;

sgtitle('Extended SEIR Model: COVID-19 with Vaccination and Mortality', ...
        'FontSize', 14, 'FontWeight', 'bold');

% --- Figure 2: Overview of all compartments together ---
figure('Position', [100, 100, 900, 450]);
plot(t, S, 'b-', 'LineWidth', 1.5); hold on;
plot(t, E, 'Color', [1 0.5 0], 'LineWidth', 1.5);
plot(t, I1, 'r-', 'LineWidth', 1.5);
plot(t, I2, 'm-', 'LineWidth', 1.5);
plot(t, R, 'g-', 'LineWidth', 1.5);
plot(t, Q, 'Color', [0.5 0 0.5], 'LineWidth', 1.5);
plot(t, H, 'Color', [0 0.5 0.5], 'LineWidth', 1.5);
plot(t, V, 'Color', [0 0.7 0], 'LineWidth', 1.5, 'LineStyle', '--');
plot(t, D, 'k-', 'LineWidth', 2);
xlabel('Time (days)'); ylabel('Number of Individuals');
title('Extended SEIR Model: All Compartments');
legend('S', 'E', 'I_1', 'I_2', 'R', 'Q', 'H', 'V', 'D', ...
       'Location', 'eastoutside');
grid on;

%% ===================== ODE FUNCTION =====================

function dydt = extended_seir(~, y, p)
    % Unpack state variables
    S  = y(1);  E  = y(2);  I1 = y(3);  I2 = y(4);
    R  = y(5);  Q  = y(6);  H  = y(7);  V  = y(8);  % D = y(9)
    
    % Total living population (excludes deceased)
    N = S + E + I1 + I2 + R + Q + H + V;
    
    % Avoid division by zero
    if N <= 0
        dydt = zeros(9,1);
        return;
    end
    
    % Force of infection
    foi = (S/N) * (p.beta1*I1 + p.beta2*I2 + p.nu*E);
    
    % ODEs
    dS  = -foi - p.q2*S + p.q1*Q + p.alpha*R - p.rho*S;
    dE  = foi - p.theta1*E - p.theta2*E;
    dI1 = p.theta1*E - p.gamma1*I1 - p.mu1*I1;
    dI2 = p.theta2*E - p.gamma2*I2 - p.varphi*I2 + p.lambda_q*(p.K + Q) - p.mu2*I2;
    dR  = p.gamma1*I1 + p.gamma2*I2 + p.phi*H - p.alpha*R;
    dH  = p.varphi*I2 - p.phi*H - p.mu3*H;
    dQ  = p.K + p.q2*S - p.lambda_q*(p.K + Q) - p.q1*Q;
    dV  = p.rho*S;
    dD  = p.mu1*I1 + p.mu2*I2 + p.mu3*H;
    
    dydt = [dS; dE; dI1; dI2; dR; dQ; dH; dV; dD];
end