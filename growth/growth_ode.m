function ydot = growth_ode(t,y,theta,L)
    
	% ODE for growth factor example
	
    RT = theta(1); k1 = theta(2); km1 = theta(3); kdeg = theta(4); kdegs = theta(5);

    ydot(1,1) = RT*kdeg - k1*L*y(1) + km1*y(2) - kdeg*y(1);
    ydot(2,1) = k1*L*y(1) - km1*y(2) - kdegs*y(2);
    
end