clear; close all; clc;


%% Population
N = 1000;


%% Parameters
% Transmission rates
beta1 = 0.5;         
beta2 = 0.15;        

% Progression rates from E
delta1 = 1/5;        % Incubation (5 days)
delta2 = 1/8;        
q_E = 0.03;          

% Recovery rates
gamma1 = 1/7;        
gamma2 = 1/10;      
gamma_Q = 1/14;      
gamma_H = 1/12;      

% Death rates
mu1 = 0.005;         
mu2 = 0.02;          
mu_Q = 0.002;        
mu_H = 0.03;        

% Hospitalisation
h = 0.2;             

% Waning immunity
omega = 1/182;       
omega_V = 1/182;     

% Vaccination
v = 0.001; 


%% Initial conditions
E0  = 5;
I1_0 = 3;
I2_0 = 1;
Q0  = 0;
H0  = 0;
R0_init = 0;
V0  = 0;
D0  = 0;
S0  = N - E0 - I1_0 - I2_0 - Q0 - H0 - R0_init - V0 - D0;

% Time
t_span = [0 365];
t_eval = linspace(0, 400, 2000);


%% Solve the ODE System
y0 = [S0; E0; I1_0; I2_0; Q0; H0; R0_init; V0; D0];

[t, y] = ode45(@(t,y) ode_extended_SEIR(t, y, beta1, beta2, ...
    delta1, delta2, q_E, v, omega, omega_V, ...
    gamma1, gamma2, gamma_Q, gamma_H, ...
    mu1, mu2, mu_Q, mu_H, h), t_eval, y0);

% Compartments
S  = y(:,1);
E  = y(:,2);
I1 = y(:,3);
I2 = y(:,4);
Q  = y(:,5);
H  = y(:,6);
R  = y(:,7);
V  = y(:,8);
D  = y(:,9);

%% Find key metrics
[peak_exposed, peak_E_idx] = max(E);
peak_E_time = t(peak_E_idx);

[peak_I1, peak_I1_idx] = max(I1);
peak_I1_time = t(peak_I1_idx);

[peak_I2, peak_I2_idx] = max(I2);
peak_I2_time = t(peak_I2_idx);

final_recovered = R(end);
final_vaccinated = V(end);
total_deaths = D(end);


%% Create figure
figure('Position', [100 100 900 600]);

hold on;
plot(t, S,  'b-',  'LineWidth', 2.5, 'DisplayName', ...
    'Susceptible');
plot(t, E,  'Color', [1 0.5 0], 'LineWidth', 2.5, ...
    'DisplayName', 'Exposed');
plot(t, I1, 'r-',  'LineWidth', 2.5, 'DisplayName', ...
    'I_1');
plot(t, I2, 'm-',  'LineWidth', 2.5, 'DisplayName', ...
    'I_2');
plot(t, R,  'g-',  'LineWidth', 2.5, 'DisplayName', ...
    'Recovered');
plot(t, Q,  'c--', 'LineWidth', 1.5, 'DisplayName', ...
    'Quarantine');
plot(t, H,  'Color', [0.5 0 0.5], 'LineStyle', '--', ...
    'LineWidth', 1.5, 'DisplayName', 'Hospitalised');
plot(t, V,  'k--', 'LineWidth', 1.5, 'DisplayName', ...
    'Vaccinated');
plot(t, D,  'Color', [0.3 0.3 0.3], 'LineStyle', ':', ...
    'LineWidth', 2, 'DisplayName', 'Deceased');


% Formatting
xlabel('Time (days)', 'FontSize', 14, 'FontWeight', ...
    'bold');
ylabel('Number of Individuals', 'FontSize', 14, ...
    'FontWeight', 'bold');
title('Extended SEIR Model: All Compartments', ...
    'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'northeast', 'FontSize', 9);
grid on;
set(gca, 'FontSize', 12);
xlim([0 400]);
ylim([0 N]);

hold off;

%% Save figure
print('fig_extended_SEIR', '-dpng', '-r300');
fprintf(['\nFigure saved as: ' ...
    '    fig_extended_SEIR.png\n']);

%% ================================================================
%  ODE function — equation
function dydt = ode_extended_SEIR(~, y, beta1, ...
    beta2,delta1, delta2, q_E, v, omega, ...
    omega_V,gamma1, gamma2, gamma_Q, ...
    gamma_H,mu1, mu2, mu_Q, mu_H, h)

    S  = y(1); 
    E  = y(2); 
    I1 = y(3); 
    I2 = y(4);
    Q  = y(5); 
    H  = y(6); 
    R  = y(7); 
    V  = y(8); 
    D  = y(9);

    N_alive = S + E + I1 + I2 + Q + H + R + V;

    dS  = -beta1*S*I1/N_alive - beta2*S*I2/N_alive - v*S + omega*R + omega_V*V;
    dE  =  beta1*S*I1/N_alive + beta2*S*I2/N_alive - (delta1 + delta2 + q_E)*E;
    dI1 =  delta1*E - (gamma1 + mu1)*I1;
    dI2 =  delta2*E - (gamma2 + mu2 + h)*I2;
    dQ  =  q_E*E - (gamma_Q + mu_Q)*Q;
    dH  =  h*I2 - (gamma_H + mu_H)*H;
    dR  =  gamma1*I1 + gamma2*I2 + gamma_Q*Q + gamma_H*H - omega*R;
    dV  =  v*S - omega_V*V;
    dD  =  mu1*I1 + mu2*I2 + mu_Q*Q + mu_H*H;

    dydt = [dS; dE; dI1; dI2; dQ; dH; dR; dV; dD];
end