function dydt = ode_SEIR(t, y, beta, sigma, gamma, N)


    S = y(1);  
    E = y(2);  
    I = y(3);  
    R = y(4); 
    
    
    % Force of Infection
    lambda = beta * I / N;
    
    
    dS = -lambda * S;              
    dE = lambda * S - sigma * E;   
    dI = sigma * E - gamma * I;    
    dR = gamma * I;
    


    
    % Return derivatives as column vector
    dydt = [dS; dE; dI; dR; dH; dQ; dD];
end