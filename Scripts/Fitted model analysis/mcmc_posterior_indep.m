% Import posterior parameter distributions for the independent transmission
% and symptoms model obtained using data augmentation MCMC, and calculate
% the posterior distribution of the proportion of presymptomatic
% transmissions

clear all; close all; clc;

addpath('../../Data')
addpath('../../Results')
addpath('../../Functions/Indep')

% Load assumed parameters

load('../../Data/assumed_parameters.mat','inc_mu','inc_sigma')

% Load output of MCMC fitting procedure

load('../../Results/param_fit_indep.mat','theta_mat','ll_vec','output')

% Calculate posterior distributions of individual model parameters

mean_post = (theta_mat(:,1));
sd_post = (theta_mat(:,2));
beta_post = (theta_mat(:,3));

% Point estimates of model parameters

mean_best = mean(mean_post);
sd_best = mean(sd_post);
beta_best = mean(beta_post);

% Posterior and point estimates for the proportion of presymptomatic
% transmissions

F_inc = @(t)logncdf(t,inc_mu,inc_sigma);

logn_mu = @(m,s)log(m.^2./sqrt(s.^2+m.^2));
logn_sigma = @(m,s)sqrt(log(1+s.^2./m.^2));

mu_post = logn_mu(mean_post,sd_post);
sigma_post = logn_sigma(mean_post,sd_post);
prob_presymp_post = get_presymp_trans_probs_indep_logn(mu_post,sigma_post,F_inc);

mu_best = logn_mu(mean_best,sd_best);
sigma_best = logn_sigma(mean_best,sd_best);
prob_presymp_best = get_presymp_trans_probs_indep_logn(mu_best,sigma_best,F_inc);

% Estimated mean and standard deviation of realised household generation
% times, and the proportion of presymptomatic transmissions

empirical_summary_mat = output.empirical_summary_mat;

% Save results

save('../../Results/mcmc_posterior_indep.mat','mean_post','mean_best','sd_post','sd_best','beta_post','beta_best','prob_presymp_post','prob_presymp_best','mu_post','sigma_post','empirical_summary_mat')

rmpath('../../Data')
rmpath('../../Results')
rmpath('../../Functions/Indep')