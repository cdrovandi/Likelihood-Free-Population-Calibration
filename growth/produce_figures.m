%% FIGURES FOR GROWTH CASE STUDY

addpath("../shared")

%% TWO UNKNOWN PARAMETERS

load('data.mat');
load('results_bsl_growth_subset.mat');

prior.num_params = 4;

prior.p1 = [2.5e5 0.25 0 0];
prior.p2 = [8e5 3 200 2];
prior.sampler = @() [unifrnd(prior.p1,prior.p2)]; 
prior.pdf = @(theta_trans) prod(exp(theta_trans)./(1 + exp(theta_trans)).^2);
prior.trans_f = @(theta) [log((theta - prior.p1)./(prior.p2 - theta))];
prior.trans_finv = @(theta_trans) [(prior.p2.*exp(theta_trans) + prior.p1)./(1 + exp(theta_trans))];

sim_func = @growth_sim_subset;
sim_params.m = 100;

thin = 20;
theta_thin = theta(1:thin:end,:);
theta_true = [6.5e5 1.7 sqrt(0.6e4) sqrt(0.05)];

M = 20000;
m = 24;

%% 2PARAM: PLOT MARGINAL POSTERIORS

figure(1); clf;
subaxis(1,4,1,'Spacing',0,'SpacingVert',0,'SpacingHoriz',0.02,'Padding',0.01,'marginL',0.03,'marginR',0.03,'mt',0.07,'mb',0.17)

subaxis(1);
[f,xi] = ksdensity(theta_thin(:,1));
plot(xi,f,'k','LineWidth',2);
hold on;
line([prior.p1(1) prior.p2(1)], [1/(prior.p2(1) - prior.p1(1)) 1/(prior.p2(1) - prior.p1(1))], 'Color', 'b', 'LineWidth',2, 'LineStyle', '--');
plot(theta_true(1), 0, 'rx','MarkerSize',10,'LineWidth',2);
xlim([6.3e5 6.7e5])
xlabel('$\mu_{R_T}$','FontSize',16)

subaxis(2);
[f,xi] = ksdensity(theta_thin(:,2));
plot(xi,f,'k','LineWidth',2);
hold on;
line([prior.p1(2) prior.p2(2)], [1/(prior.p2(2) - prior.p1(2)) 1/(prior.p2(2) - prior.p1(2))], 'Color', 'b', 'LineWidth',2, 'LineStyle', '--');
plot(theta_true(2), 0, 'rx','MarkerSize',10,'LineWidth',2);
xlim([1.4 1.9])
xlabel('$\mu_{k_1}$','FontSize',16)

subaxis(3);
[f,xi] = ksdensity(theta_thin(:,3));
plot(xi,f,'k','LineWidth',2);
hold on;
line([prior.p1(3) prior.p2(3)], [1/(prior.p2(3) - prior.p1(3)) 1/(prior.p2(3) - prior.p1(3))], 'Color', 'b', 'LineWidth',2, 'LineStyle', '--');
plot(theta_true(3), 0, 'rx','MarkerSize',10,'LineWidth',2);
xlim([0 200])
xlabel('$\sigma_{R_T}$','FontSize',16)

subaxis(4);
[f,xi] = ksdensity(theta_thin(:,4));
plot(xi,f,'k','LineWidth',2);
hold on;
line([prior.p1(4) prior.p2(4)], [1/(prior.p2(4) - prior.p1(4)) 1/(prior.p2(4) - prior.p1(4))], 'Color', 'b', 'LineWidth',2, 'LineStyle', '--');
plot(theta_true(4), 0, 'rx','MarkerSize',10,'LineWidth',2);
xlim([0 0.5])
xlabel('$\sigma_{k_1}$','FontSize',16)


cleanfigure
matlab2tikz('growth_posterior_theta.tex','parseStrings',false,...
    'width','\textwidth',...
    'height','0.23\textwidth',...
    'extraCode','\pgfplotsset{tick label style={font=\footnotesize}}')


%% 2PARAM: POSTERIOR PREDICTIVE SUMMARIES (COMPUTE)

[M,~] = size(theta_thin);
ssx = zeros(M,5);

sim_func = @growth_sim_subset;
sim_params.m = 100;
parfor i = 1:M
    x = sim_func(theta_thin(i,:), sim_params);
    covx = cov(x);
    ssx(i,:) = [mean(x) std(x) covx(1,2)];
end

covy = cov(y);
ssy = [mean(y) std(y) covy(1,2)];


%% 2PARAM: PLOT POSTERIOR PREDICTIVE SUMMARIES

figure(2); clf;
subaxis(1,5,1,'Spacing',0,'SpacingVert',0,'SpacingHoriz',0.02,'Padding',0.01,'marginL',0.03,'marginR',0.03,'mt',0.07,'mb',0.14)

subaxis(1);
[f,xi] = ksdensity(ssx(:,1));
plot(xi,f,'k','LineWidth',2);
hold on;
plot(ssy(1), 0, 'rx','MarkerSize',10,'LineWidth',2);
xlabel('$S_1$','FontSize',16)

subaxis(2);
[f,xi] = ksdensity(ssx(:,2));
plot(xi,f,'k','LineWidth',2);
hold on;
plot(ssy(2), 0, 'rx','MarkerSize',10,'LineWidth',2);
xlabel('$S_2$','FontSize',16)

subaxis(3);
[f,xi] = ksdensity(ssx(:,3));
plot(xi,f,'k','LineWidth',2);
hold on;
plot(ssy(3), 0, 'rx','MarkerSize',10,'LineWidth',2);
xlabel('$S_3$','FontSize',16);
xticks([1000,1400,1800])
xlim([800,1800])

subaxis(4);
[f,xi] = ksdensity(ssx(:,4));
plot(xi,f,'k','LineWidth',2);
hold on;
plot(ssy(4), 0, 'rx','MarkerSize',10,'LineWidth',2);
xlabel('$S_4$','FontSize',16)

subaxis(5);
[f,xi] = ksdensity(ssx(:,5));
plot(xi,f,'k','LineWidth',2);
hold on;
plot(ssy(5), 0, 'rx','MarkerSize',10,'LineWidth',2);
xlabel('$S_5$','FontSize',16)

cleanfigure
matlab2tikz('growth_posterior_output.tex','parseStrings',false,...
    'width','\textwidth',...
    'height','0.2\textwidth',...
    'extraCode','\pgfplotsset{tick label style={font=\footnotesize},x tick label style={/pgf/number format/precision=5}}')

%% 2PARAM: PLOT POSTERIORS OF f(x)
figure(3); clf;
subaxis(1,2,1,'Spacing',0,'SpacingVert',0,'SpacingHoriz',0.02,'Padding',0.01,'marginL',0.03,'marginR',0.03,'mt',0.07,'mb',0.14)

subaxis(1)

min1 = theta_true(1) - 5*theta_true(3);
max1 = theta_true(1) + 5*theta_true(3);
%grid1 = min1:((max1-min1)/100):max1;
grid1 = linspace(6.496e5,6.506e5,101);

hold on;
for i = 1:M
    out1(i,:) = normpdf(grid1,theta_thin(i,1),theta_thin(i,3));
end
fill_between(grid1,quantile(out1,0.025),quantile(out1,0.975))
plot(grid1,median(out1),'Color',[0.6,0.6,0.6],'LineWidth',2);
plot(grid1,normpdf(grid1,theta_true(1),theta_true(3)),'r--','LineWidth',2);

ylim([0,0.007]);
xlabel('$R_T$');
legend('95\% PI','Median','True');
box on;

subaxis(2)
min2 = theta_true(2) - 5*theta_true(4);
max2 = theta_true(2) + 5*theta_true(4);
%grid2 = min2:((max2-min2)/100):max2;
grid2 = linspace(0.6,3.0,101);

hold on;
for i = 1:M
    out2(i,:) = normpdf(grid2,theta_thin(i,2),theta_thin(i,4));
end
fill_between(grid2,quantile(out2,0.025),quantile(out2,0.975))
plot(grid2,median(out2),'Color',[0.6,0.6,0.6],'LineWidth',2);
plot(grid2,normpdf(grid2,theta_true(2),theta_true(4)),'r--','LineWidth',2);

xlabel('$k_1$');
legend('95\% PI','Median','True');
box on;

cleanfigure
matlab2tikz('growth_posterior_input.tex','parseStrings',false,...
    'width','\textwidth',...
    'height','0.3\textwidth',...
    'extraCode','\pgfplotsset{legend style={font=\footnotesize},tick label style={font=\footnotesize}}')

%% 2PARAM: PLOT POSTERIORS SAMPLES FOR R_T DISTRIBUTION
figure(4); clf;

grid1 = linspace(theta_true(1)-4*theta_true(3),theta_true(1)+4*theta_true(3),101);

theta_thin2 = theta_thin(50:50:end,:);

hold on; 
plot(grid1,normpdf(grid1,theta_true(1),theta_true(3)),'r--','LineWidth',2);
for i = 1:size(theta_thin2,1)
    m = theta_thin2(i,1);
    s = theta_thin2(i,3);
    xgrid = linspace(m-4*s,m+4*s,101);
    plot(xgrid, normpdf(xgrid,m,s), 'Color', [0.5 0.5 0.5]);
end
plot(grid1,normpdf(grid1,theta_true(1),theta_true(3)),'r--','LineWidth',2);
ylim([0,0.02])
xlabel('$R_T$');
legend("True","Posterior sample");
box on;

matlab2tikz('growth_posterior_sims_input.tex','parseStrings',false,...
    'width','0.985\textwidth',...
    'height','0.2\textwidth',...
    'extraCode','\pgfplotsset{tick label style={font=\footnotesize}}')


%% FIVE UNKNOWN PARAMETERS
load('results_bsl_growth.mat');

sim_func = @growth_sim;
sim_params.m = 100;
prior.p1 = [2.5e5 0.25 2 0.005 0.1 0 0 0 0 0];
prior.p2 = [8e5 3 20 0.1 0.5 200 2 2 2 2];

thin = 50;
theta_thin = theta(1:thin:end,:);
[M,~] = size(theta_thin);
theta_true = [6.5e5 1.7 8 0.015 0.25 77.4597 0.2236 0.0001 0.0001 0.0001];

%% 5PARAM: MCMC CHAINS
figure(5); clf;

subaxis(2,5,1,'Spacing',0,'SpacingVert',0,'SpacingHoriz',0.02,'Padding',0.01,'marginL',0.03,'marginR',0.03,'mt',0.03,'mb',0,'PaddingBottom',0.08);

labels = {'$\mu_{R_T}$', '$\mu_{k_1}$', '$\mu_{k_{-1}}$', '$\mu_{k_\text{deg}}$', '$\mu_{k_\text{deg}^*}$'...,
    '$\sigma_{R_T}$', '$\sigma_{k_1}$', '$\sigma_{k_{-1}}$', '$\sigma_{k_\text{deg}}$', '$\sigma_{k_\text{deg}^*}$'};

for i = 1:10
    subaxis(i);
    plot(theta_thin(:,i),'k')
    hold on;
    plot(1:M, repmat(theta_true(i),1,M), 'r', 'LineWidth',2);
    xticks([0,500,1000])
    xtickangle(0)
    title(labels{i})
end


cleanfigure
matlab2tikz('growth_theta_trace_5param.tex','parseStrings',false,...
    'width','\textwidth',...
    'height','0.4\textwidth',...
    'extraCode','\pgfplotsset{legend style={font=\scriptsize},scaled x ticks=false, x tick label style={/pgf/number format/fixed}, tick label style={font=\scriptsize}}')

%% 5PARAM: POSTERIOR PREDICTIVE SUMMARIES (COMPUTE)

[M,~] = size(theta_thin);
ssx = zeros(M,5);

sim_func = @growth_sim;
sim_params.m = 100;
parfor i = 1:M
    x = sim_func(theta_thin(i,:), sim_params);
    covx = cov(x);
    ssx(i,:) = [mean(x) std(x) covx(1,2)];
end

covy = cov(y);
ssy = [mean(y) std(y) covy(1,2)];


%% 5PARAM: PLOT POSTERIOR PREDICTIVE SUMMARIES

figure(6);
subaxis(1,5,1,'Spacing',0,'SpacingVert',0,'SpacingHoriz',0.02,'Padding',0.01,'marginL',0.03,'marginR',0.03,'mt',0.07,'mb',0.14)

subaxis(1);
[f,xi] = ksdensity(ssx(:,1));
plot(xi,f,'k','LineWidth',2);
hold on;
plot(ssy(1), 0, 'rx','MarkerSize',10,'LineWidth',2);
xlabel('$S_1$','FontSize',16)

subaxis(2);
[f,xi] = ksdensity(ssx(:,2));
plot(xi,f,'k','LineWidth',2);
hold on;
plot(ssy(2), 0, 'rx','MarkerSize',10,'LineWidth',2);
xlabel('$S_2$','FontSize',16)

subaxis(3);
[f,xi] = ksdensity(ssx(:,3));
plot(xi,f,'k','LineWidth',2);
hold on;
plot(ssy(3), 0, 'rx','MarkerSize',10,'LineWidth',2);
xlabel('$S_3$','FontSize',16)
xticks([1000,1400,1800])
xlim([800,1800])
xtickangle(0)

subaxis(4);
[f,xi] = ksdensity(ssx(:,4));
plot(xi,f,'k','LineWidth',2);
hold on;
plot(ssy(4), 0, 'rx','MarkerSize',10,'LineWidth',2);
xlabel('$S_4$','FontSize',16)

subaxis(5);
[f,xi] = ksdensity(ssx(:,5));
plot(xi,f,'k','LineWidth',2);
hold on;
plot(ssy(5), 0, 'rx','MarkerSize',10,'LineWidth',2);
xlabel('$S_5$','FontSize',16)

cleanfigure
matlab2tikz('growth_posterior_output_5param.tex','parseStrings',false,...
    'width','\textwidth',...
    'height','0.2\textwidth',...
    'extraCode','\pgfplotsset{tick label style={font=\footnotesize}}')

%% 5PARAM: POSTERIOR PREDICTIVE DISTRIBUTIONS
figure(7); clf;
labels = {'$R_T$', '$k_1$', '$k_{-1}$', '$k_\text{deg}$', '$k_\text{deg}^*$'};
r = 10;
ind = randsample(1:M, r, 'false');

figure;
subaxis(1,5,1,'Spacing',0,'SpacingVert',0,'SpacingHoriz',0.02,'Padding',0.01,'marginL',0.03,'marginR',0.03,'mt',0.04,'mb',0.16)

for i = 1:5
    subaxis(i);
    hold on;
    for j = 1:r
        [f,xi] = ksdensity(normrnd(theta(ind(j),i), theta(ind(j),i+5), 10000, 1));
        plot(xi,f,'k')
    end
    xlabel(labels{i}, 'FontSize', 16);
end

cleanfigure
matlab2tikz('growth_input_sample_5param.tex','parseStrings',false,...
    'width','\textwidth',...
    'height','0.2\textwidth',...
    'extraCode','\pgfplotsset{tick label style={font=\footnotesize}}')