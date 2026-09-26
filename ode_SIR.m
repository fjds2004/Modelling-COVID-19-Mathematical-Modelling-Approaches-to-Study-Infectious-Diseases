function dydt = ode_SIR(t, y, beta, gamma, N)

    S = y(1);  % Susceptible
    I = y(2);  % Infected
    R = y(3);  % Recovered
    
    % Calculate force of infection
    lambda = beta * I / N;
    
    % System of differential equations
    dS = -lambda * S;           % Susceptibles becoming infected
    dI = lambda * S - gamma * I; % New infections minus recoveries
    dR = gamma * I;              % Recoveries
    
    % Return derivatives as column vector
    dydt = [dS; dI; dR];
end