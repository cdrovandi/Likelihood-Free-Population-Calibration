



load('data.mat');

% define prior
prior.num_params = 5;
prior.sampler = @() [sort(normrnd(0,1,1,2)) exprnd(1) exprnd(1) rand]; 
prior.pdf = @(theta_trans) prior_pdf_twocomp(theta_trans);
prior.trans_f = @(theta) [theta(1) theta(2) log(theta(3)) log(theta(4)) log(theta(5)/(1-theta(5)))];
prior.trans_finv = @(theta_trans) [theta_trans(1) theta_trans(2) exp(theta_trans(3)) exp(theta_trans(4)) 1/(1+exp(-theta_trans(5)))];

% define model
sim_func = @normal_twocomp;
sim_params.m = 1000;


M = 100000; % number of MCMC iterations
m = 50; % number of model simulations used to estimate synthetic likelhiood
start = prior.trans_f([0.3 0.5 0.015 0.043 1/3]);

% obtained from pilot runs
cov_rw = [3.78606543920221e-06,1.71217912938229e-06,0.000619560467825367,-5.43705557027783e-05,5.09892571190340e-05;1.71217912938229e-06,2.09301717928309e-06,0.000457514688261072,-3.39906076622041e-05,3.48107268118504e-05;0.000619560467825367,0.000457514688261072,0.540105237638016,-0.00918086901540768,0.0105785818174206;-5.43705557027783e-05,-3.39906076622041e-05,-0.00918086901540768,0.00187328075757633,-0.00111159398602610;5.09892571190340e-05,3.48107268118504e-05,0.0105785818174206,-0.00111159398602610,0.00145413841327850];
numComp=2;

[theta, loglike] = bayes_bsl_aux(y,m,M,start,cov_rw,prior,sim_func,sim_params,numComp);
save('Results_BSL.mat','theta','loglike');



