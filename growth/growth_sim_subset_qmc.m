function h = growth_sim_subset_qmc(theta,sim_params)

% function to generate population data in the growth factor example
% in this function the parameters k_-1, k_deg and k_deg^* are fixed

mu = theta(1:2); sigma = theta(3:4);
m = sim_params.m;
r1 = sim_params.r1;
r2 = sim_params.r2;
h = zeros(m,2);

parfor i = 1:m
    thetam = norminv(r1(i,:),mu,sigma);
    [~,sol1] = ode45(@(t,y)growth_ode(t,y,[thetam 8 0.015 0.25],2),[0 10],[0;0]);
    %h(i,1) = sol(end,2);
    thetam = norminv(r2(i,:),mu,sigma);
    [~,sol2] = ode45(@(t,y)growth_ode(t,y,[thetam 8 0.015 0.25],10),[0 10],[0;0]);
    %h(i,2) = sol(end,2);
    h(i,:) = [sol1(end,2) sol2(end,2)];
end

end