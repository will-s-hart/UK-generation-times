% Produce the panels in Figure 1-figure supplement 2 of our manuscript.

% This script the export_fig package (freely available at
% https://github.com/altmany/export_fig) to export the figure panels to pdf
% files.

clear all; close all; clc;

addpath('../Data')
addpath('../Results')

% Load results

load('../Results/mcmc_posterior_mech.mat','p_E_post','mu_inv_post','alpha_post','beta_post')

% Colorscheme

c1 = [0, 0.4470, 0.7410]; c2 = [0.8500, 0.3250, 0.0980];
c3 = [0.9290, 0.6940, 0.1250]; c4 = [0.4940, 0.1840, 0.5560];

% Pre-format figures

for k = 1:4
figsetup(k)
end

% Ratio between the shape parameters of the Gamma distributed latent (E)
% and incubation (combined E and P) periods, p_E=k_E/k_inc

figure(1); hold on;
p2 = histogram(p_E_post,'facecolor',c1,'normalization','pdf','edgecolor','auto');
xlim([0,1])
p1 = plot(0:0.001:1,betapdf(0:0.001:1,2.1,2.1),'k--','linewidth',3);
p3 = plot(mean(p_E_post)*[1,1],[min(ylim),max(ylim)],'r--','linewidth',3);
xlabel({'{\itk}_{E}/{\itk}_{inc}'})
ylabel('Density')
l = legend([p1,p2,p3],{'Prior','Posterior','Posterior mean'});
l.FontSize = 15;
l.Position = [0.5160    0.8040    0.2860    0.1180];

% Mean symptomatic infectious (I) period, mu_inv=1/mu

figure(2); hold on;
p2 = histogram(mu_inv_post,'facecolor',c1,'normalization','pdf','edgecolor','auto');
xlim([0,10])
p1 = plot(0:0.01:10,lognpdf(0:0.01:10,1.6,0.8),'k--','linewidth',3);
p3 = plot(mean(mu_inv_post)*[1,1],[min(ylim),max(ylim)],'r--','linewidth',3);
xlabel({'1/\mu (days)'})
ylabel('Density')
l = legend([p1,p2,p3],{'Prior','Posterior','Posterior mean'});
l.FontSize = 15;
l.Position = [0.5160    0.8040    0.2860    0.1180];

% Ratio of the transmission rates during the presymptomatic infectious (P)
% and symptomatic infectious (I) periods, alpha=alpha_P

figure(3); hold on;
p2 = histogram(alpha_post,'facecolor',c1,'normalization','pdf','edgecolor','auto');
xlim([0,8])
p1 = plot(0:0.01:8,lognpdf(0:0.01:8,0,0.8),'k--','linewidth',3);
p3 = plot(mean(alpha_post)*[1,1],[min(ylim),max(ylim)],'r--','linewidth',3);
ylim([0,0.7])
yticks(0:0.1:0.7)
xlabel({'\alpha_P'})
ylabel('Density')
l = legend([p1,p2,p3],{'Prior','Posterior','Posterior mean'});
l.FontSize = 15;
l.Position = [0.5160    0.8040    0.2860    0.1180];

% Overall infectiousness parameter, beta_0

figure(4); hold on;
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

for k = 1:4
figsetup(k)
end

% Export figure panels to pdf

figure(1); export_fig Figures/figure1_supp2/figure1_supp2A.pdf -nocrop -painters -transparent
figure(2); export_fig Figures/figure1_supp2/figure1_supp2B.pdf -nocrop -painters -transparent
figure(3); export_fig Figures/figure1_supp2/figure1_supp2C.pdf -nocrop -painters -transparent
figure(4); export_fig Figures/figure1_supp2/figure1_supp2D.pdf -nocrop -painters -transparent

rmpath('../Data')
rmpath('../Results')