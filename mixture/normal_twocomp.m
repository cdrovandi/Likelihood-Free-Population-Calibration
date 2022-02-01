function h = normal_twocomp(theta,sim_params)

	% function to simulate population data

    mu1 = theta(1); mu2 = theta(2); sigma1 = theta(3); sigma2 = theta(4); w = theta(5); 
    m = sim_params.m;

    x = zeros(m,1);
    ind = rand(m,1)<w;
    
    x(ind) = normrnd(mu1, sigma1, sum(ind==1),1);
    x(~ind) = normrnd(mu2, sigma2, sum(ind==0),1);

    h = normrnd(x, 0.045);

end