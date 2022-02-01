%% FIGURES FOR INTERNALISATION/FLOW CYTOMETRY CASE STUDY
addpath("../shared")

%% LOAD JLD2 FILE FROM GITHUB AND CONVERT TO .MAT
%  vv run this to obtain the .mat file using the "results_to_matlab.jl"
%  script.
%!/Applications/Julia-1.7.app/Contents/Resources/julia/bin/julia "results_to_matlab.jl"
load("Browning2021.mat")

%% PLOT MARGINAL POSTERIORS

thin = 100;
theta_thin = theta(1:thin:end,:);

param_names_latex = ["$\alpha_1$","$\alpha_2$","$\sigma_1$","$\sigma_2$","$\mu_R$","$\sigma_R$","$\mu_\lambda$","$\sigma_\lambda$","$\omega_\lambda$","$\mu_\beta$","$\sigma_\beta$","$\omega_\beta$","$\rho_{R\lambda}$","$\rho_{R\beta}$","$\bar{\rho}_{\lambda\beta}$","$p$"];

figure(1); clf();
subaxis(3,6,1,'Spacing',0,'SpacingVert',0.1,'SpacingHoriz',0.02,'Padding',0.01,'marginL',0.03,'marginR',0.03,'mt',0.05,'mb',0.13)

for i = 1:16
    subaxis(i);
    [f,xi] = ksdensity(theta_thin(:,i));
    plot(xi,f,'k','LineWidth',2);
    xlabel(param_names_latex(i),'FontSize',16);
end

cleanfigure;
matlab2tikz('internalisation_posterior.tex','parseStrings',false,...
    'width','\textwidth',...
    'height','0.7\textwidth',...
    'extraCode','\pgfplotsset{legend style={font=\footnotesize}, tick label style={font=\footnotesize}, scaled y ticks=true}')

%% PLOT POSTERIORS OF DISTRIBUTIONS

figure(2); clf;
subaxis(1,3,1,'Spacing',0,'SpacingVert',0.1,'SpacingHoriz',0.02,'Padding',0.01,'marginL',0.03,'marginR',0.03,'mt',0.05,'mb',0.13)

x1 = linspace(0,3,200);
x2 = linspace(0,0.5,200);
x3 = linspace(0,0.25,200);
f1 = zeros(size(theta_thin,1),length(x1));
f2 = zeros(size(theta_thin,1),length(x2));
f3 = zeros(size(theta_thin,1),length(x3));
for i = 1:size(theta_thin,1)
    f1(i,:) = altlognormpdf(x1,theta_thin(i,5),theta_thin(i,6));
    f2(i,:) = altgampdf(x2,theta_thin(i,7),theta_thin(i,8),theta_thin(i,9));
    f3(i,:) = altgampdf(x3,theta_thin(i,10),theta_thin(i,11),theta_thin(i,12));
end

subaxis(1); hold on; box on;
fill_between(x1,quantile(f1,0.025),quantile(f1,0.975));
plot(x1,median(f1),'Color',[0.6,0.6,0.6],'LineWidth',2);
plot(x1,altlognormpdf(x1,theta_best(5),theta_best(6)),'b--','LineWidth',2);
xlabel("$R$"); legend("95\% CI","Median","Best fit")

subaxis(2); hold on; box on;
fill_between(x2,quantile(f2,0.025),quantile(f2,0.975));
plot(x2,median(f2),'Color',[0.6,0.6,0.6],'LineWidth',2);
plot(x2,altgampdf(x2,theta_best(7),theta_best(8),theta_best(9)),'b--','LineWidth',2);
xlabel("$\lambda$"); legend("95\% CI","Median","Best fit")

subaxis(3); hold on; box on;
fill_between(x3,quantile(f3,0.025),quantile(f3,0.975));
plot(x3,median(f3),'Color',[0.6,0.6,0.6],'LineWidth',2);
plot(x3,altgampdf(x3,theta_best(10),theta_best(11),theta_best(12)),'b--','LineWidth',2);
xlabel("$\beta$"); legend("95\% CI","Median","Best fit")

cleanfigure;
matlab2tikz('internalisation_predictive.tex','parseStrings',false,...
    'width','\textwidth',...
    'height','0.3\textwidth',...
    'extraCode','\pgfplotsset{legend style={font=\footnotesize}, tick label style={font=\footnotesize}}')

%% POSTERIOR PREDICTIVE DISTRIBUTION

%  vv run this to obtain the .mat file using the "InternalisationPredictive.jl"
%  script.
%!/Applications/Julia-1.7.app/Contents/Resources/julia/bin/julia "InternalisationPredictive.jl"
load("Browning2021_Predictive.mat")
figure(3); clf;
subaxis(1,4,1,'Spacing',0,'SpacingVert',0.1,'SpacingHoriz',0.02,'Padding',0.01,'marginL',0.03,'marginR',0.03,'mt',0.05,'mb',0.13)

subaxis(1); hold on; box on;
fill_between(x',quantile(h1',0.025),quantile(h1',0.975))
plot(x,median(h1'),'Color',[0.6,0.6,0.6],'LineWidth',2);
plot(x,B(:,1),'b--','LineWidth',2)
ylim([0,10])
xlabel("Proportion")
title("1 min")

subaxis(2); hold on; box on;
fill_between(x',quantile(h2',0.025),quantile(h2',0.975))
plot(x,median(h2'),'Color',[0.6,0.6,0.6],'LineWidth',2);
plot(x,B(:,2),'b--','LineWidth',2)
ylim([0,5])
xlabel("Proportion")
title("5 min")

subaxis(3); hold on; box on;
fill_between(x',quantile(h3',0.025),quantile(h3',0.975))
plot(x,median(h3'),'Color',[0.6,0.6,0.6],'LineWidth',2);
plot(x,B(:,3),'b--','LineWidth',2)
ylim([0,12])
xlabel("Proportion")
title("30 min")

subaxis(4); hold on; box on;
fill_between(x',quantile(h4',0.025),quantile(h4',0.975))
plot(x,median(h4'),'Color',[0.6,0.6,0.6],'LineWidth',2);
plot(x,B(:,4),'b--','LineWidth',2)
xlim([0.9,1.0])
ylim([0,100])
xlabel("Proportion")
title("180 min")

cleanfigure;
matlab2tikz('internalisation_predictive_proportion.tex','parseStrings',false,...
    'width','\textwidth',...
    'height','0.3\textwidth',...
    'extraCode','\pgfplotsset{legend style={font=\footnotesize}, tick label style={font=\footnotesize}}')
