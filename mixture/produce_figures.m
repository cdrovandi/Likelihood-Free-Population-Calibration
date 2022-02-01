%% FIGURES FOR GROWTH CASE STUDY

addpath("../shared")

%% LOAD RESULTS

load('data.mat');
load('Results_BSL.mat');

thin = 100;
M = size(theta,1);

sim_func = @normal_twocomp;
sim_params.m = 1000;

theta_thin = theta(1:thin:end,:);
N = M/thin;

%% PLOT MARGINAL POSTERIORS

theta_true = [0.3 0.5 0.015 0.043 1/3];

figure(1); clf;
subaxis(1,5,1,'Spacing',0,'SpacingVert',0,'SpacingHoriz',0.02,'Padding',0.01,'marginL',0.03,'marginR',0.03,'mt',0.05,'mb',0.13)

subaxis(1);
[f,xi] = ksdensity(theta_thin(:,1));
plot(xi,f,'k','LineWidth',2);
hold on;
plot(theta_true(1), 0, 'rx','MarkerSize',10,'LineWidth',2);
xlabel('\mu_1','FontSize',16);

subaxis(2);
[f,xi] = ksdensity(theta_thin(:,2));
plot(xi,f,'k','LineWidth',2);
hold on;
plot(theta_true(2), 0, 'rx','MarkerSize',10,'LineWidth',2);
xlabel('\mu_2','FontSize',16);

subaxis(3);
[f,xi] = ksdensity(theta_thin(:,3));
plot(xi,f,'k','LineWidth',2);
hold on;
plot(theta_true(3), 0, 'rx','MarkerSize',10,'LineWidth',2);
xlabel('\sigma_1','FontSize',16);

subaxis(4);
[f,xi] = ksdensity(theta_thin(:,4));
plot(xi,f,'k','LineWidth',2);
hold on;
plot(theta_true(4), 0, 'rx','MarkerSize',10,'LineWidth',2);
xlabel('\sigma_2','FontSize',16);

subaxis(5);
[f,xi] = ksdensity(theta_thin(:,5));
plot(xi,f,'k','LineWidth',2);
hold on;
plot(theta_true(5), 0, 'rx','MarkerSize',10,'LineWidth',2);
xlabel('\omega','FontSize',16);

cleanfigure
matlab2tikz('mixture_posteriors.tex',...
    'width','\textwidth',...
    'height','0.2\textwidth',...
    'extraCode','\pgfplotsset{scaled x ticks=false, x tick label style={/pgf/number format/fixed}, tick label style={font=\footnotesize}}')


%% PLOT POSTERIOR PREDICTIVE DISTRIBUTION

xgrid = 0:0.01:1;
f = zeros(N,length(xgrid));
figure(2); clf;
hold on;
for i = 1:N
   ys = sim_func(theta_thin(i,:), sim_params);
   [f(i,:),~] = ksdensity(ys,xgrid);
end
[fobs,~] = ksdensity(y,xgrid);
fill_between(xgrid,quantile(f,0.025),quantile(f,0.975))
plot(xgrid,median(f),'Color',[0.6,0.6,0.6],'LineWidth',2);
plot(xgrid,fobs,'r--','LineWidth',2);
xlim([0 1])
ylim([0,5])
legend('95\% PI','Median','Data')
set(gca,'XTick',0:0.1:1);
set(gca,'YTick',0:1:5);
xlabel("$y$")
ylabel("$h(y)$")
box on;

cleanfigure
matlab2tikz('mixture_predictive.tex','parseStrings',false,...
    'width','0.4\textwidth',...
    'height','0.25\textwidth',...
    'extraCode','\pgfplotsset{legend style={font=\footnotesize},scaled x ticks=false, x tick label style={/pgf/number format/fixed}, tick label style={font=\footnotesize}}')

%% PLOT POSTERIOR OF f(x)

xgrid = linspace(0.2,0.65,500);
f = zeros(N,length(xgrid));
figure(2); clf;
hold on;
f = zeros(N,length(xgrid));
for i = 1:(M/thin)
    f(i,:) = theta_thin(i,5)*normpdf(xgrid, theta_thin(i,1), theta_thin(i,3)) + (1 - theta_thin(i,5))*normpdf(xgrid, theta_thin(i,2), theta_thin(i,4));
end
[fobs,~] = ksdensity(y,xgrid);
fill_between(xgrid,quantile(f,0.025),quantile(f,0.975))
plot(xgrid,median(f),'Color',[0.6,0.6,0.6],'LineWidth',2);
plot(xgrid, 1/3*normpdf(xgrid, 0.3, 0.015) + (1 - 1/3)*normpdf(xgrid, 0.5, 0.043),'r--','LineWidth',2);
xlim([0.2 0.65])
ylim([0.0,15.0])
legend('95\% PI','Median','True')
set(gca,'XTick',0:0.1:1);
set(gca,'YTick',0:5:15);
xlabel("$x$")
ylabel("$f(x)$")
box on;

cleanfigure
matlab2tikz('mixture_posteriorfx.tex','parseStrings',false,...
    'width','0.4\textwidth',...
    'height','0.25\textwidth',...
    'extraCode','\pgfplotsset{legend style={font=\footnotesize},scaled x ticks=false, x tick label style={/pgf/number format/fixed}, tick label style={font=\footnotesize}}')

