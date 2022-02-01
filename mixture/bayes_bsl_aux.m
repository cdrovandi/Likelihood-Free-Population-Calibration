function [theta, loglike] = bayes_bsl_aux(y,m,M,start,cov_rw,prior,sim_func,sim_params,numComp)
%%
% function to run BSL on mixture example
%
% inputs:
% y - population data
% m - number of model simulations for estimating synthetic likelihood
% M - number of MCMC iterations
% start - initial parameter values of chain on transformed space
% cov_rw - MCMC random walk covariance on transformed space
% prior - contains functions related to the prior distribution
% sim_func - function to simulate population data from the model
% sim_params - holds auxiliary objects that may be required by sim_func
% numComp - number of components in the mixture used to summarise the population data
%
% outputs:
% theta - MCMC chain of parameter values on transformed space
% loglike - MCMC chain of synthetic log-likelihood values
%%


obj = gmdistribution.fit(y,numComp,'Replicates',100,'Options',statset('MaxIter',10000,'TolFun',1e-10));

theta_d = [obj.PComponents(1:(numComp-1)) obj.mu' reshape(obj.Sigma,numComp,1)'];
ssy = zeros(1,length(theta_d));

ssx = zeros(m,length(theta_d));
for k = 1:m
    x = sim_func(prior.trans_finv(start), sim_params);
    ssx(k,:) = compute_grad(theta_d,x,obj,numComp);
end


%Calculating the mean and covariance of the summary statistics
the_mean = mean(ssx);
the_cov = cov(ssx);
L = chol(the_cov);
logdetA = 2*sum(log(diag(L)));

% synthetic likelihood
loglike_curr = -0.5*logdetA - 0.5*(ssy-the_mean)*inv(the_cov)*(ssy-the_mean)';


theta = zeros(M,prior.num_params);
loglike = zeros(M,1);

% MH - IL
theta_curr = start;

for i = 1:M
    
    theta_prop = mvnrnd(theta_curr,cov_rw);
    
    logprior_curr = log(prior.pdf(theta_curr));
    logprior_prop = log(prior.pdf(theta_prop));
    
    prop = prior.trans_finv(theta_prop);
    
    ssx = zeros(m,length(theta_d));
    for k = 1:m
        x = sim_func(prop, sim_params);
        ssx(k,:) = compute_grad(theta_d,x,obj,numComp);
    end
    if (any(any(isnan(ssx))))
        theta(i,:) = theta_curr;
        loglike(i) = loglike_curr;
        continue;
    end
        
    
    %Calculating the mean and covariance of the summary statistics
    the_mean = mean(ssx);
    the_cov = cov(ssx);
    [L,p] = chol(the_cov);
    if (p>0)
        theta(i,:) = theta_curr;
        loglike(i) = loglike_curr;
        continue;
    end
        
    logdetA = 2*sum(log(diag(L)));
    
    % synthetic likelihood
    loglike_prop = -0.5*logdetA - 0.5*(ssy-the_mean)*inv(the_cov)*(ssy-the_mean)';

    mh = exp(loglike_prop - loglike_curr + logprior_prop - logprior_curr);
    
    if (mh > rand)
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
