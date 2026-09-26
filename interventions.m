clear; clc; close all;

%% UK second wave parameters
N = 66739900;       % UK population
T_days = 250;       % approx second wave 

% --- Epidemiological parameters
% Transmission rates
beta1 = 0.35;       
beta2 = 0.12;       

% Progression rates from E
delta1 = 1/5.2;     % incubation (5.2 days)
delta2 = 1/8.0;     % slower, identified individuals

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
h_rate = 0.3;       

% Waning immunity
omega = 1/180;      % 6 months
omega_V = 1/180;    

% --- Control parameters ---
v_vax = 0.003;      % approx 200k/day at peak
q_E_val = 0.05;     % quarantine rate

%% Initial conditions (fractions of N)
I1_0 = 3000 / N;
E_0  = 2 * I1_0;
I2_0 = I1_0 * 0.15;
H_0  = 500 / N;
Q_0  = 0;
R_0  = 0.05;        % 5% immune from first wave
V_0  = 0;
D_0  = 41000 / N;   % cumulative deaths from first wave
S_0  = 1 - E_0 - I1_0 - I2_0 - Q_0 - H_0 - R_0 - V_0 - D_0;

y0 = [S_0, E_0, I1_0, I2_0, Q_0, H_0, R_0, V_0, D_0];

tspan = [0 T_days];

%% Run three scenarios
% No intervention (v=0, q_E=0)
params1 = [beta1, beta2, delta1, delta2, 0, 0, omega, omega_V, ...
           gamma1, gamma2, gamma_Q, gamma_H, mu1, mu2, mu_Q, mu_H, h_rate];
[t1, y1] = ode45(@(t,y) extended_model(t, y, params1), tspan, y0);

% Vaccination only (v>0, q_E=0)
params2 = [beta1, beta2, delta1, delta2, 0, v_vax, omega, omega_V, ...
           gamma1, gamma2, gamma_Q, gamma_H, mu1, mu2, mu_Q, mu_H, h_rate];
[t2, y2] = ode45(@(t,y) extended_model(t, y, params2), tspan, y0);

% Combined (v>0, q_E>0)
params3 = [beta1, beta2, delta1, delta2, q_E_val, v_vax, omega, omega_V, ...
           gamma1, gamma2, gamma_Q, gamma_H, mu1, mu2, mu_Q, mu_H, h_rate];
[t3, y3] = ode45(@(t,y) extended_model(t, y, params3), tspan, y0);

%% Figure: side-by-side plots
figure('Position', [100 100 1200 500]);

% (a) Infected population I1 + I2
subplot(1,2,1);
plot(t1, (y1(:,3) + y1(:,4)) * N, 'r-', 'LineWidth', 2); hold on;
plot(t2, (y2(:,3) + y2(:,4)) * N, 'Color', [1 0.5 0], ...
     'LineStyle', '--', 'LineWidth', 2);
plot(t3, (y3(:,3) + y3(:,4)) * N, 'g-.', 'LineWidth', 2);
xlabel('Time (Days)', 'FontSize', 12);
ylabel('Active Infections (I_1 + I_2)', 'FontSize', 12);
title('Infected Population Under Different Scenarios', 'FontSize', 13);
legend('No Intervention', 'Vaccination Only', 'Combined (Vax + Q)', ...
       'Location', 'northwest', 'FontSize', 10);
grid on; grid minor;

% (b) Total deaths
subplot(1,2,2);
plot(t1, (y1(:,9) - y1(1,9)) * N, 'r-', 'LineWidth', 2); hold on;
plot(t2, (y2(:,9) - y2(1,9)) * N, 'Color', [1 0.5 0], ...
     'LineStyle', '--', 'LineWidth', 2);
plot(t3, (y3(:,9) - y3(1,9)) * N, 'g-.', 'LineWidth', 2);
xlabel('Time (Days)', 'FontSize', 12);
ylabel('Cumulative Deaths (Second Wave)', 'FontSize', 12);
title('Total Deaths Under Different Scenarios', 'FontSize', 13);
legend('No Intervention', 'Vaccination Only', 'Combined (Vax + Q)', ...
       'Location', 'northwest', 'FontSize', 10);
grid on; grid minor;

%% Extended model ODE function 
function dydt = extended_model(~, y, p)
    S = y(1); E = y(2); I1 = y(3); I2 = y(4);
    Q = y(5); H = y(6); R = y(7); V = y(8); D = y(9);

    b1 = p(1);  b2 = p(2);
    d1 = p(3);  d2 = p(4);  qE = p(5);
    v  = p(6);  om = p(7);  omV = p(8);
    g1 = p(9);  g2 = p(10); gQ = p(11); gH = p(12);
    m1 = p(13); m2 = p(14); mQ = p(15); mH = p(16);
    hr = p(17);

    N_alive = S + E + I1 + I2 + Q + H + R + V;

    dS  = -b1*S*I1/N_alive - b2*S*I2/N_alive - v*S + om*R + omV*V;
    dE  =  b1*S*I1/N_alive + b2*S*I2/N_alive - (d1 + d2 + qE)*E;
    dI1 =  d1*E - (g1 + m1)*I1;
    dI2 =  d2*E - (g2 + m2 + hr)*I2;
    dQ  =  qE*E - (gQ + mQ)*Q;
    dH  =  hr*I2 - (gH + mH)*H;
    dR  =  g1*I1 + g2*I2 + gQ*Q + gH*H - om*R;
    dV  =  v*S - omV*V;
    dD  =  m1*I1 + m2*I2 + mQ*Q + mH*H;

    dydt = [dS; dE; dI1; dI2; dQ; dH; dR; dV; dD];
end