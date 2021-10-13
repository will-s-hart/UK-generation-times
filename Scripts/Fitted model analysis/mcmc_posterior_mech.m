% Import posterior parameter distributions for the mechanistic model
% obtained using data augmentation MCMC, and calculate the posterior
% distributions of the mean and standard deviation of generation times and
% the proportion of presymptomatic transmissions

clear all; close all; clc;

addpath('../../Data')
addpath('../../Results')
addpath('../../Functions/Mech')

% Load assumed parameters

load('../../Data/assumed_parameters.mat','k_inc','gamma','k_I','params_known')

% Load output of MCMC fitting procedure

load('../../Results/param_fit_mech.mat','theta_mat','ll_vec','output')

% Calculate posterior distributions of individual model parameters

p_E_post = (theta_mat(:,1));
k_E_post = k_inc*p_E_post;
k_P_post = k_inc - k_E_post;
mu_inv_post = (theta_mat(:,2));
alpha_post = (theta_mat(:,3));
beta_post = (theta_mat(:,4));

% Point estimates of model parameters

p_E_best = mean(p_E_post);
k_E_best = k_inc*p_E_best;
k_P_best = k_inc - k_E_best;
mu_inv_best = mean(mu_inv_post);
alpha_best = mean(alpha_post);
beta_best = mean(beta_post);

theta_best = [p_E_best,mu_inv_best,(alpha_best),(beta_best)]; %fitted parameters

params_best = get_params_mech(theta_best,params_known); %all parameters

% Posterior and point estimates for the proportion of presymptomatic
% transmissions

prob_presymp_post = (alpha_post.*k_P_post/(k_inc*gamma))./((alpha_post.*k_P_post/(k_inc*gamma))+mu_inv_post);
prob_presymp_best = (alpha_best*k_P_best/(k_inc*gamma))/((alpha_best*k_P_best/(k_inc*gamma))+mu_inv_best);

% Posterior and point estimates for the mean and standard deviation of
% generation times

no_steps_kept = length(k_E_post);
k_inc_post = repmat(params_known(1),no_steps_kept,1);
gamma_post = repmat(params_known(2),no_steps_kept,1);
k_I_post = repmat(params_known(3),no_steps_kept,1);
rho_post = repmat(params_known(4),no_steps_kept,1);
x_A_post = repmat(params_known(5),no_steps_kept,1);
params_post = [gamma_post,1./mu_inv_post,k_inc_post,k_E_post,k_I_post,alpha_post,beta_post,rho_post,x_A_post];

[mean_post,sd_post] = get_gen_mean_sd_mech(params_post);
[mean_best,sd_best] = get_gen_mean_sd_mech(params_best);

% Save results

empirical_summary_mat = output.empirical_summary_mat;

save('../../Results/mcmc_posterior_mech.mat','p_E_post','mu_inv_post','alpha_post','beta_post','prob_presymp_post','params_best','mean_post','mean_best','sd_post','sd_best','empirical_summary_mat')

rmpath('../../Data')
rmpath('../../Results')
rmpath('../../Functions/Mech')