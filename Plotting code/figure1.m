% Produce the panels in Figure 1 of our manuscript.

% This script requires the Violinplot-Matlab package to produce the violin
% plots (freely available at https://github.com/bastibe/Violinplot-Matlab),
% and requires the export_fig package (freely available at
% https://github.com/altmany/export_fig) to export the figure panels to pdf
% files.

clear all; close all; clc;

addpath('../Data')
addpath('../Results')

% Load results

load('../Results/mcmc_posterior_indep.mat','mean_post','sd_post','beta_post','prob_presymp_post')

mean_indep = mean_post; sd_indep = sd_post;
beta_indep = beta_post; presymp_indep = prob_presymp_post;

load('../Results/mcmc_posterior_mech.mat','mean_post','sd_post','beta_post','prob_presymp_post')

mean_mech = mean_post; sd_mech = sd_post;
beta_mech = beta_post; presymp_mech = prob_presymp_post;

% Colorscheme

c1 = [0, 0.4470, 0.7410]; c2 = [0.8500, 0.3250, 0.0980];
c3 = [0.9290, 0.6940, 0.1250]; c4 = [0.4940, 0.1840, 0.5560];

% Pre-format figures

for k = 1:4
figsetup(k)
end

% Mean generation time

figure(1); hold on;
v = violinplot([mean_indep,mean_mech],{'',''},'ShowData',false,'ViolinAlpha',1);
v(1).ViolinColor = c1; v(2).ViolinColor = c2;
v1 = v(1).ViolinPlot; v2 = v(2).ViolinPlot;
ylim([0,10])
ylabel({'Mean generation time (days)'})
l = legend([v1,v2],{'Independent transmission/symptoms','Mechanistic'});
l.FontSize = 15;
l.Position = [0.2200    0.8070    0.5800    0.0810];

% Standard deviation of generation times

figure(2); hold on;
v = violinplot([sd_indep,sd_mech],{'',''},'ShowData',false,'ViolinAlpha',1);
v(1).ViolinColor = c1; v(2).ViolinColor = c2;
v1 = v(1).ViolinPlot; v2 = v(2).ViolinPlot;
ylim([0,10])
ylabel({'Standard deviation of generation';'times (days)'})
l = legend([v1,v2],{'Independent transmission/symptoms','Mechanistic'});
l.FontSize = 15;
l.Position = [0.2200    0.8070    0.5800    0.0810];

% Proportion of transmissions before symptom onset

figure(3); hold on;
v = violinplot([presymp_indep,presymp_mech],{'',''},'ShowData',false,'ViolinAlpha',1);
v(1).ViolinColor = c1; v(2).ViolinColor = c2;
v1 = v(1).ViolinPlot; v2 = v(2).ViolinPlot;
ylim([0,1])
ylabel({'Proportion of transmissions';'before symptom onset'})
l = legend([v1,v2],{'Independent transmission/symptoms','Mechanistic'});
l.FontSize = 15;
l.Position = [0.2200    0.8070    0.5800    0.0810];

% Overall infectiousness parameter, beta_0

figure(4); hold on;
v = violinplot([beta_indep,beta_mech],{'',''},'ShowData',false,'ViolinAlpha',1);
v(1).ViolinColor = c1; v(2).ViolinColor = c2;
v1 = v(1).ViolinPlot; v2 = v(2).ViolinPlot;
ylim([0,2.5])
ylabel({'Overall infectiousness, \beta_0'})
l = legend([v1,v2],{'Independent transmission/symptoms','Mechanistic'});
l.FontSize = 15;
l.Position = [0.2200    0.8070    0.5800    0.0810];

% Post-format figures

for k = 1:4
figsetup(k)
end

% Export figure panels to pdf

figure(1); export_fig Figures/figure1/figure1A.pdf -nocrop -painters -transparent
figure(2); export_fig Figures/figure1/figure1B.pdf -nocrop -painters -transparent
figure(3); export_fig Figures/figure1/figure1C.pdf -nocrop -painters -transparent
figure(4); export_fig Figures/figure1/figure1D.pdf -nocrop -painters -transparent

rmpath('../Data')
rmpath('../Results')