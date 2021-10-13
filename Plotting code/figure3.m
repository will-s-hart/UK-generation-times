% Produce the panels in Figure 3 of our manuscript.

% This script requires the Violinplot-Matlab package to produce the violin
% plots (freely available at https://github.com/bastibe/Violinplot-Matlab),
% and requires the export_fig package (freely available at
% https://github.com/altmany/export_fig) to export the figure panels to pdf
% files.

clear all; close all; clc;

addpath('../Data')
addpath('../Results')

% Load results

load('../Results/Split_data/mcmc_posterior_indep_splitdata.mat','mean_post1','mean_post2','mean_post3','sd_post1','sd_post2','sd_post3','beta_post1','beta_post2','beta_post3','prob_presymp_post1','prob_presymp_post2','prob_presymp_post3')

% Colorscheme

c1 = [0, 0.4470, 0.7410]; c2 = [0.8500, 0.3250, 0.0980];
c3 = [0.9290, 0.6940, 0.1250]; c4 = [0.4940, 0.1840, 0.5560];

% Pre-format figures

for k = 1:4
figsetup(k)
end

% Mean generation time

figure(1); hold on;
v = violinplot([mean_post1,mean_post2,mean_post3],{'','',''},'ShowData',false,'ViolinAlpha',1);
v(1).ViolinColor = c1; v(2).ViolinColor = c2; v(3).ViolinColor = c3;
v1 = v(1).ViolinPlot; v2 = v(2).ViolinPlot; v3 = v(3).ViolinPlot;
ylim([0,10])
ylabel({'Mean generation time (days)'})
l = legend([v1,v2,v3],{'March-April','May-August','September-November'},'location','n');
l.FontSize = 15;
l.Position = [0.4170    0.8060    0.3820    0.1180];

% Standard deviation of generation times

figure(2); hold on;
inds = (sd_post1<=15&sd_post2<=15&sd_post3<=15); %remove outliers to ensure smooth violin
v = violinplot([sd_post1(inds),sd_post2(inds),sd_post3(inds)],{'','',''},'ShowData',false,'ViolinAlpha',1);
v(1).ViolinColor = c1; v(2).ViolinColor = c2; v(3).ViolinColor = c3;
v1 = v(1).ViolinPlot; v2 = v(2).ViolinPlot; v3 = v(3).ViolinPlot;
ylim([0,10])
ylabel({'Standard deviation of generation';'times (days)'})
l = legend([v1,v2,v3],{'March-April','May-August','September-November'},'location','n');
l.FontSize = 15;
l.Position = [0.4170    0.8060    0.3820    0.1180];

% Proportion of transmissions before symptom onset

figure(3); hold on;
v = violinplot([prob_presymp_post1,prob_presymp_post2,prob_presymp_post3],{'','',''},'ShowData',false,'ViolinAlpha',1);
v(1).ViolinColor = c1; v(2).ViolinColor = c2; v(3).ViolinColor = c3;
v1 = v(1).ViolinPlot; v2 = v(2).ViolinPlot; v3 = v(3).ViolinPlot;
ylim([0,1])
ylabel({'Proportion of transmissions';'before symptom onset'})
l = legend([v1,v2,v3],{'March-April','May-August','September-November'},'location','n');
l.FontSize = 15;
l.Position = [0.4170    0.8060    0.3820    0.1180];

% Overall infectiousness parameter, beta_0

figure(4); hold on;
v = violinplot([beta_post1,beta_post2,beta_post3],{'','',''},'ShowData',false,'ViolinAlpha',1);
v(1).ViolinColor = c1; v(2).ViolinColor = c2; v(3).ViolinColor = c3;
v1 = v(1).ViolinPlot; v2 = v(2).ViolinPlot; v3 = v(3).ViolinPlot;
ylim([0,3])
ylabel({'Overall infectiousness, \beta_0'})
l = legend([v1,v2,v3],{'March-April','May-August','September-November'},'location','n');
l.FontSize = 15;
l.Position = [0.4170    0.8060    0.3820    0.1180];

% Post-format figures

for k = 1:4
figsetup(k)
end

% Export figure panels to pdf

figure(1); export_fig Figures/figure3/figure3A.pdf -nocrop -painters -transparent
figure(2); export_fig Figures/figure3/figure3B.pdf -nocrop -painters -transparent
figure(3); export_fig Figures/figure3/figure3C.pdf -nocrop -painters -transparent
figure(4); export_fig Figures/figure3/figure3D.pdf -nocrop -painters -transparent

rmpath('../Data')
rmpath('../Results')