clear; clc; close all;


%% UK second wave parameters
N = 66739900;
T_days = 250;

%% Initial parameters
baseline = [0.35, % beta1
            0.12, % beta2
            1/5.2, % delta1
            1/8.0, % delta2
            0.05, % q_E
            0.003, % v
            1/180, % omega
            1/180, % omega_V
            1/7, % gamma1
            1/10, % gamma2
            1/14, % gamma_Q
            1/12, % gamma_H
            0.005, % mu1
            0.02, % mu2
            0.002, % mu_Q
            0.03, % mu_H
            0.3];% h

% Varying Parameters
param_names = {'\beta_1', '\delta_1', '\gamma_1', 'v', 'q_E', 'h'};
param_indices = [1, 3, 9, 6, 5, 17]; 

perturbation = 0.10;  % +/- 10%


%% Initial conditions (fractions of N)
I1_0 = 3000 / N;
E_0  = 2 * I1_0;
I2_0 = I1_0 * 0.15;
H_0  = 500 / N;
Q_0  = 0;
R_0  = 0.05;
V_0  = 0;
D_0  = 41000 / N;
S_0  = 1 - E_0 - I1_0 - I2_0 - Q_0 - H_0 - R_0 - V_0 - D_0;

y0 = [S_0, E_0, I1_0, I2_0, Q_0, H_0, R_0, V_0, D_0];
tspan = [0 T_days];

%% Base values
baseline_cum = compute_cumulative_infections(baseline, y0, tspan, N);


%% OAT Sensitivity Analysis
n_params = length(param_indices);
low_pct  = zeros(n_params, 1);
high_pct = zeros(n_params, 1);

fprintf('--- Results ---\n');
for k = 1:n_params
    idx = param_indices(k);

    % Low perturbation (-10%)
    p_low = baseline;
    p_low(idx) = baseline(idx) * (1 - perturbation);
    cum_low = compute_cumulative_infections(p_low, y0, tspan, N);

    % High perturbation (+10%)
    p_high = baseline;
    p_high(idx) = baseline(idx) * (1 + perturbation);
    cum_high = compute_cumulative_infections(p_high, y0, tspan, N);

    low_pct(k)  = (cum_low - baseline_cum) / baseline_cum * 100;
    high_pct(k) = (cum_high - baseline_cum) / baseline_cum * 100;

    fprintf('%8s: -10%% -> %12.0f (%+.1f%%)  | +10%% -> %12.0f (%+.1f%%)\n', ...
            param_names{k}, cum_low, low_pct(k), cum_high, high_pct(k));
end

%% Tornado Plot
% Sort by total swing
total_swing = abs(high_pct - low_pct);
[~, sort_idx] = sort(total_swing, 'ascend');

sorted_names = param_names(sort_idx);
sorted_low   = low_pct(sort_idx);
sorted_high  = high_pct(sort_idx);

figure('Position', [100 100 1000 600]);
y_pos = 1:n_params;

for k = 1:n_params
    %  -10% perturbation
    barh(y_pos(k), sorted_low(k), 0.55, 'FaceColor', [1 0 0], ...
         'FaceAlpha', 0.85, 'EdgeColor', 'w');
    hold on;
    % +10% perturbation
    barh(y_pos(k), sorted_high(k), 0.55, 'FaceColor', [0 0 1], ...
         'FaceAlpha', 0.85, 'EdgeColor', 'w');
end

% Formatting
set(gca, 'YTick', y_pos, 'YTickLabel', sorted_names, 'FontSize', 13);
xlabel('Percentage Change in Total Infections (%)', 'FontSize', 12);
title('Sensitivity of Total Infections to \pm10% Parameter Variation', ...
      'FontSize', 13, 'FontWeight', 'bold');
xline(0, 'k-', 'LineWidth', 1);
grid on;


%% Helper function: compute cumulative infections

function cum_inf = compute_cumulative_infections(params, y0, tspan, N)
    [t, y] = ode45(@(t,y) extended_model(t, y, params), tspan, y0);

    % Cumulative infections = integral of (delta1 + delta2 + q_E) * E(t)
    d1 = params(3);
    d2 = params(4);
    qE = params(5);
    outflow_E = (d1 + d2 + qE) .* y(:,2); 
    cum_inf = trapz(t, outflow_E) * N;
end

%% Extended model ODE function — equation (1)

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