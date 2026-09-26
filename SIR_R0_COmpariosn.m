clear; close all; clc;


N = 1000;           
gamma = 0.2;        
I0 = 1;             
R0_init = 0;        
S0 = N - I0;        

t_span = [0 200];
t_eval = linspace(100, 200, 1000);

%% Different R0 values
R0_values = [0.8, 1.5, 2.5, 4.0];

colors = [0 0 1;   
          1 0 0;   
          0 1 0;   
          0 0 0];  

fprintf('Comparing SIR Dynamics for Different R0 Values\n');
fprintf('================================================\n\n');

% Store results
results = struct();

figure('Position', [100 100 800 600]);
hold on;

for i = 1:length(R0_values)
    R0_val = R0_values(i);
    beta = R0_val * gamma;
    
    y0 = [S0; I0; R0_init];
    [t, y] = ode45(@(t,y) ode_SIR(t, y, beta, gamma, N), t_eval, y0);
    
    I = y(:,2);
    [peak, peak_idx] = max(I);
    peak_time = t(peak_idx);
    
    % Store results
    results(i).R0 = R0_val;

    results(i).peak = peak;

    results(i).peak_time = peak_time;

    results(i).final_size = y(end,3);

    results(i).t = t;

    results(i).I = I;
    

    %% Plot
    if R0_val < 1
        plot(t, I, '--', 'LineWidth', 2.5, 'Color', colors(i,:), ...
             'DisplayName', sprintf('R_0 = %.1f (no epidemic)', R0_val));
    else
        plot(t, I, '-', 'LineWidth', 2.5, 'Color', colors(i,:), ...
             'DisplayName', sprintf('R_0 = %.1f', R0_val));
    end
end

xlabel('Time (days)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Infected Individuals', 'FontSize', 14, 'FontWeight', 'bold');
title('Effect of R_0 on Epidemic Peak', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'northeast', 'FontSize', 10);
grid on;
set(gca, 'FontSize', 11);
hold off;

fprintf('\nFigure saved as: fig_SIR_endemic_R0_varying.png\n');