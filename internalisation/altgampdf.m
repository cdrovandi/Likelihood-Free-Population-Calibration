function density = altgampdf(x,mu,sigma,omega)

    % Normal parameterisation
    alpha = 4 / omega^2;
    theta = sigma * abs(omega) / 2;
    
    % Shift and negative scale
    if omega > 0
        shift = alpha * theta - mu;
        d = truncate(makedist('Gamma','a',alpha,'b',theta),shift,inf);
    else
        shift = alpha * theta + mu;
        d = truncate(makedist('Gamma','a',alpha,'b',theta),-inf,shift);
    end
    
    density = pdf(d,sign(omega) * x + shift);

end