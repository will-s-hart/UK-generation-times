% Import posterior parameter distributions obtained by fitting the
% independent transmission and symptoms model to data from three
% sub-intervals of the study period using data augmentation MCMC, and
% calculate the posterior distribution of the proportion of presymptomatic
% transmissions

clear all; close all; clc;

addpath('../../Data')
addpath('../../Results')
addpath('../../Functions/Indep')

% Load assumed parameters

load('../../Data/assumed_parameters.mat','inc_mu','inc_sigma')

% Load output of MCMC fitting procedure

load('../../Results/Split_data/param_fit_indep_splitdata.mat','theta_mat1','ll_vec1','output1','theta_mat2','ll_vec2','output2','theta_mat3','ll_vec3','output3')

% Calculate posterior distributions of individual model parameters

mean_post1 = (theta_mat1(:,1));
sd_post1 = (theta_mat1(:,2));
beta_post1 = (theta_mat1(:,3));

mean_post2 = (theta_mat2(:,1));
sd_post2 = (theta_mat2(:,2));
beta_post2 = (theta_mat2(:,3));

mean_post3 = (theta_mat3(:,1));
sd_post3 = (theta_mat3(:,2));
beta_post3 = (theta_mat3(:,3));

% Posterior proportion of presymptomatic transmissions

F_inc = @(t)logncdf(t,inc_mu,inc_sigma);

logn_mu = @(m,s)log(m.^2./sqrt(s.^2+m.^2));
logn_sigma = @(m,s)sqrt(log(1+s.^2./m.^2));

mu_post1 = logn_mu(mean_post1,sd_post1);
sigma_post1 = logn_sigma(mean_post1,sd_post1);
prob_presymp_post1 = get_presymp_trans_probs_indep_logn(mu_post1,sigma_post1,F_inc);

mu_post2 = logn_mu(mean_post2,sd_post2);
sigma_post2 = logn_sigma(mean_post2,sd_post2);
prob_presymp_post2 = get_presymp_trans_probs_indep_logn(mu_post2,sigma_post2,F_inc);

mu_post3 = logn_mu(mean_post3,sd_post3);
sigma_post3 = logn_sigma(mean_post3,sd_post3);
prob_presymp_post3 = get_presymp_trans_probs_indep_logn(mu_post3,sigma_post3,F_inc);

% Save results

save('../../Results/Split_data/mcmc_posterior_indep_splitdata.mat','mean_post1','mean_post2','mean_post3','sd_post1','sd_post2','sd_post3','beta_post1','beta_post2','beta_post3','prob_presymp_post1','prob_presymp_post2','prob_presymp_post3')

rmpath('../../Data')
rmpath('../../Results')
rmpath('../../Functions/Indep')