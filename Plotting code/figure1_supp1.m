% Produce the panels in Figure 1-figure supplement 1 of our manuscript.

% This script the export_fig package (freely available at
% https://github.com/altmany/export_fig) to export the figure panels to pdf
% files.

clear all; close all; clc;

addpath('../Data')
addpath('../Results')

% Load results

load('../Results/mcmc_posterior_indep.mat','mean_post','sd_post','beta_post')

% Colorscheme

c1 = [0, 0.4470, 0.7410]; c2 = [0.8500, 0.3250, 0.0980];
c3 = [0.9290, 0.6940, 0.1250]; c4 = [0.4940, 0.1840, 0.5560];

% Pre-format figures

for k = 1:3
figsetup(k)
end

% Mean generation time

figure(1); hold on;
p2 = histogram(mean_post,'facecolor',c1,'normalization','pdf','edgecolor','auto');
xlim([0,10])
p1 = plot(0:0.01:10,lognpdf(0:0.01:10,1.6,0.35),'k--','linewidth',3);
p3 = plot(mean(mean_post)*[1,1],[min(ylim),max(ylim)],'r--','linewidth',3);
xlabel({'Mean generation time (days)'})
ylabel('Density')
l = legend([p1,p2,p3],{'Prior','Posterior','Posterior mean'});
l.FontSize = 15;
l.Position = [0.5160    0.8040    0.2860    0.1180];

% Standard deviation of generation times

figure(2); hold on;
p2 = histogram(sd_post,'facecolor',c1,'normalization','pdf','edgecolor','auto');
xlim([0,10])
p1 = plot(0:0.01:10,lognpdf(0:0.01:10,0.7,0.65),'k--','linewidth',3);
p3 = plot(mean(sd_post)*[1,1],[min(ylim),max(ylim)],'r--','linewidth',3);
xlabel({'Standard deviation of generation';'times (days)'})
ylabel('Density')
l = legend([p1,p2,p3],{'Prior','Posterior','Posterior mean'});
l.FontSize = 15;
l.Position = [0.5160    0.8040    0.2860    0.1180];

% Overall infectiousness parameter, beta_0

figure(3); hold on;
p2 = histogram(beta_post,'facecolor',c1,'normalization','pdf','edgecolor','auto');
xlim([0,2.5])
p1 = plot(0:0.01:2.5,lognpdf(0:0.01:2.5,0.7,0.8),'k--','linewidth',3);
p3 = plot(mean(beta_post)*[1,1],[min(ylim),max(ylim)],'r--','linewidth',3);
xlabel({'\beta_0'})
ylabel('Density')
l = legend([p1,p2,p3],{'Prior','Posterior','Posterior mean'});
l.FontSize = 15;
l.Position = [0.5160    0.8040    0.2860    0.1180];

% Post-format figures

for k = 1:3
figsetup(k)
end

% Export figure panels to pdf

figure(1); export_fig Figures/figure1_supp1/figure1_supp1A.pdf -nocrop -painters -transparent
figure(2); export_fig Figures/figure1_supp1/figure1_supp1B.pdf -nocrop -painters -transparent
figure(3); export_fig Figures/figure1_supp1/figure1_supp1C.pdf -nocrop -painters -transparent

rmpath('../Data')
rmpath('../Results')