function h = growth_sim_subset(theta,sim_params)

	% function to generate population data in the growth factor example
	% in this function the parameters k_-1, k_deg and k_deg^* are fixed

    mu = theta(1:2); sigma = theta(3:4);
    m = sim_params.m;
    h = zeros(m,2); 
    
    for i = 1:m
       thetam = normrnd(mu,sigma);
       [~,sol] = ode45(@(t,y)growth_ode(t,y,[thetam 8 0.015 0.25],2),[0 10],[0;0]);
       h(i,1) = sol(end,2);
       thetam = normrnd(mu,sigma);
       [~,sol] = ode45(@(t,y)growth_ode(t,y,[thetam 8 0.015 0.25],10),[0 10],[0;0]);
       h(i,2) = sol(end,2);
    end

end