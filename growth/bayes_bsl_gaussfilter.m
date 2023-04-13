function [theta, loglike] = bayes_bsl_gaussfilter(y,M,start,cov_rw,prior,sim_func,sim_params)
%%
% function to run filter inference on growth factor example
%
% inputs:
% y - population data
% M - number of MCMC iterations
% start - initial parameter values of chain on transformed space
% cov_rw - MCMC random walk covariance on transformed space
% prior - contains functions related to the prior distribution
% sim_func - function to simulate population data from the model
% sim_params - holds auxiliary objects that may be required by sim_func
%
% outputs:
% theta - MCMC chain of parameter values on transformed space
% loglike - MCMC chain of synthetic log-likelihood values
%%



% run model simulations and compute moments
x = sim_func(prior.trans_finv(start), sim_params);

%Calculating the mean and covariance of data
the_mean = mean(x);
the_cov = cov(x);

% synthetic likelihood
loglike_curr = sum(log(mvnpdf(y,the_mean,the_cov)));

% store chain of values
theta = zeros(M,prior.num_params);
loglike = zeros(M,1);

theta_curr = start;

for i = 1:M
    i
	% propose new parameter value on transformed space
    theta_prop = mvnrnd(theta_curr,cov_rw);
    
    logprior_curr = log(prior.pdf(theta_curr));
    logprior_prop = log(prior.pdf(theta_prop));
    
	% back transform to simulate
    prop = prior.trans_finv(theta_prop);
    
    % run model simulations and compute moments
    x = sim_func(prop, sim_params);
    if (any(any(isnan(x))))
        theta(i,:) = theta_curr;
        loglike(i) = loglike_curr;
        continue;
    end
        
    %Calculating the mean and covariance of the summary statistics
    the_mean = mean(x);
    the_cov = cov(x);
    
    % synthetic likelihood
    loglike_prop = sum(log(mvnpdf(y,the_mean,the_cov)));

	% MH for accepting/rejecting proposal
    mh = exp(loglike_prop - loglike_curr + logprior_prop - logprior_curr);
    
    if (mh > rand)
		% accept so update values
        theta_curr = theta_prop;
        loglike_curr = loglike_prop;
    end
    theta(i,:) = theta_curr;
    loglike(i) = loglike_curr;
    
end

% back transform
for i=1:M
    theta(i,:) = prior.trans_finv(theta(i,:));
end


end
