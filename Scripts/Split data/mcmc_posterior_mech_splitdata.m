% Import posterior parameter distributions obtained by fitting the
% mechanistic model to data from three sub-intervals of the study period
% using data augmentation MCMC, and calculate the posterior distributions
% of the mean and standard deviation of generation times and the proportion
% of presymptomatic transmissions

clear all; close all; clc;

addpath('../../Data')
addpath('../../Results')
addpath('../../Functions/Mech')

% Load assumed parameters

load('../../Data/assumed_parameters.mat','k_inc','gamma','k_I','params_known')

% Load output of MCMC fitting procedure

load('../../Results/Split_data/param_fit_mech_splitdata.mat','theta_mat1','ll_vec1','output1','theta_mat2','ll_vec2','output2','theta_mat3','ll_vec3','output3')

% Calculate posterior distributions of individual model parameters

p_E_post1 = (theta_mat1(:,1));
k_E_post1 = k_inc*p_E_post1;
k_P_post1 = k_inc - k_E_post1;
mu_inv_post1 = (theta_mat1(:,2));
alpha_post1 = (theta_mat1(:,3));
beta_post1 = (theta_mat1(:,4));

p_E_post2 = (theta_mat2(:,1));
k_E_post2 = k_inc*p_E_post2;
k_P_post2 = k_inc - k_E_post2;
mu_inv_post2 = (theta_mat2(:,2));
alpha_post2 = (theta_mat2(:,3));
beta_post2 = (theta_mat2(:,4));

p_E_post3 = (theta_mat3(:,1));
k_E_post3 = k_inc*p_E_post3;
k_P_post3 = k_inc - k_E_post3;
mu_inv_post3 = (theta_mat3(:,2));
alpha_post3 = (theta_mat3(:,3));
beta_post3 = (theta_mat3(:,4));

% Posterior estimates for the mean and standard deviation of generation
% times

no_steps_kept = length(k_E_post1);
k_inc_post = repmat(params_known(1),no_steps_kept,1);
gamma_post = repmat(params_known(2),no_steps_kept,1);
k_I_post = repmat(params_known(3),no_steps_kept,1);
rho_post = repmat(params_known(4),no_steps_kept,1);
x_A_post = repmat(params_known(5),no_steps_kept,1);

params_post1 = [gamma_post,1./mu_inv_post1,k_inc_post,k_E_post1,k_I_post,alpha_post1,beta_post1,rho_post,x_A_post];
[mean_post1,sd_post1] = get_gen_mean_sd_mech(params_post1);

params_post2 = [gamma_post,1./mu_inv_post2,k_inc_post,k_E_post2,k_I_post,alpha_post2,beta_post2,rho_post,x_A_post];
[mean_post2,sd_post2] = get_gen_mean_sd_mech(params_post2);

params_post3 = [gamma_post,1./mu_inv_post3,k_inc_post,k_E_post3,k_I_post,alpha_post3,beta_post3,rho_post,x_A_post];
[mean_post3,sd_post3] = get_gen_mean_sd_mech(params_post3);

% Posterior estimates for the proportion of transmissions before symptom
% onset

prob_presymp_post1 = (alpha_post1.*k_P_post1/(k_inc*gamma))./((alpha_post1.*k_P_post1/(k_inc*gamma))+mu_inv_post1);
prob_presymp_post2 = (alpha_post2.*k_P_post2/(k_inc*gamma))./((alpha_post2.*k_P_post2/(k_inc*gamma))+mu_inv_post2);
prob_presymp_post3 = (alpha_post3.*k_P_post3/(k_inc*gamma))./((alpha_post3.*k_P_post3/(k_inc*gamma))+mu_inv_post3);

% Save results

save('../../Results/Split_data/mcmc_posterior_mech_splitdata.mat','mean_post1','mean_post2','mean_post3','sd_post1','sd_post2','sd_post3','beta_post1','beta_post2','beta_post3','prob_presymp_post1','prob_presymp_post2','prob_presymp_post3')

rmpath('../../Data')
rmpath('../../Results')
rmpath('../../Functions/Mech')