function h = growth_sim_qmc(theta,sim_params)

	% function to generate population data in the growth factor example
    % this function uses common QMC numbers

    mu = theta(1:5); sigma = theta(6:10);
    m = sim_params.m;
    r1 = sim_params.r1;
    r2 = sim_params.r2;
    h = zeros(m,2);
    
    parfor i = 1:m
       thetam = norminv(r1(i,:),mu,sigma);
       [~,sol1] = ode45(@(t,y)growth_ode(t,y,thetam,2),[0 10],[0;0]);
       %h(i,1) = sol(end,2);
       thetam = norminv(r2(i,:),mu,sigma);
       [~,sol2] = ode45(@(t,y)growth_ode(t,y,thetam,10),[0 10],[0;0]);
       %h(i,2) = sol(end,2);
       h(i,:) = [sol1(end,2) sol2(end,2)];
    end

end