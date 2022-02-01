function h = growth_sim(theta,sim_params)

	% function to generate population data in the growth factor example

    mu = theta(1:5); sigma = theta(6:10);
    m = sim_params.m;
    h = zeros(m,2); 
    
    for i = 1:m
       thetam = normrnd(mu,sigma);
       [~,sol] = ode45(@(t,y)growth_ode(t,y,thetam,2),[0 10],[0;0]);
       h(i,1) = sol(end,2);
       thetam = normrnd(mu,sigma);
       [~,sol] = ode45(@(t,y)growth_ode(t,y,thetam,10),[0 10],[0;0]);
       h(i,2) = sol(end,2);
    end

end