% Fit the parameters of the independent transmission and symptoms model to
% the data from three sub-intervals of the study interval (the suffixes
% 1,2,3 in variable names refer to the different time-intervals), using
% data augmentation MCMC

% The vector of unknown model parameters, theta  = [m_gen,s_gen,beta_0], is
% estimated in the model fitting procedure. Here, m_gen and s_gen are the
% mean and standard deviation of the generation time distribution,
% respectively, and beta_0 is the overall infectiousness parameter

clear all; close all; clc;

addpath('../../Data')
addpath('../../Functions/MCMC')
addpath('../../Functions/Indep')

% Initialise random number generator

rng(4)

% Load data

load('../../Data/data_split.mat','data_split')

data_struct_observed1 = data_split{1};
data_struct_observed2 = data_split{2};
data_struct_observed3 = data_split{3};

% Load incubation period distribution function and assumed parameters

load('../../Data/assumed_parameters.mat','f_inc_logn','rho','x_A')
f_inc = f_inc_logn;

% Number of steps in chain, keeping only some steps after burn-in and
% thinning

no_steps = 10000000;
steps_keep = ((no_steps/5)+1):100:no_steps;

% Function handles giving the density and cumululative distribution
% functions of the generation time for parameters theta

logn_mu = @(m,s)log(m^2/sqrt(s^2+m^2));
logn_sigma = @(m,s)sqrt(log(1+s^2/m^2));

f_gen_form = @(t_gen,theta)lognpdf(t_gen,logn_mu((theta(1)),(theta(2))),logn_sigma((theta(1)),(theta(2))));
F_gen_form = @(t_gen,theta)logncdf(t_gen,logn_mu((theta(1)),(theta(2))),logn_sigma((theta(1)),(theta(2))));

% Function handle giving the log-likelihood for parameters theta and
% augmented data data_struct_augmented

ll_household_form = @(theta,data_struct_augmented)log_likelihood_household_indep(f_inc,(theta(3)),rho,x_A,@(t_g)f_gen_form(t_g,theta),@(t_g)F_gen_form(t_g,theta),data_struct_augmented);

% Function handle giving the estimated mean and standard deviation of
% realised household generation times, and the presymptomatic proportion of
% transmissions, for parameters theta and augmented data
% data_struct_augmented

empirical_summary_form = @(theta,data_struct_augmented)empirical_summary_indep((theta(3)),rho,x_A,@(t_g)f_gen_form(t_g,theta),data_struct_augmented);

% Standard deviations of proposal distributions for individual model
% parameters

sd_prop_mean1 = 0.425;
sd_prop_sd1 = 2*sd_prop_mean1;
sd_prop_beta1 = 0.3*sd_prop_mean1;

sd_prop_mean2 = 0.65;
sd_prop_sd2 = 1.7*sd_prop_mean2;
sd_prop_beta2 = 0.2*sd_prop_mean2;

sd_prop_mean3 = 0.325;
sd_prop_sd3 = 3.5*sd_prop_mean3;
sd_prop_beta3 = 0.35*sd_prop_mean3;

% Correlation between proposal distributions

corr_prop_mat1 = eye(3);
corr_prop_mat2 = eye(3);
corr_prop_mat3 = eye(3);

% Covariance matrix of multivariate normal proposal distribution for the
% vector of fitted model parameters, theta

sd_prop_diag1 = diag([sd_prop_mean1,sd_prop_sd1,sd_prop_beta1]);
theta_prop_cov_mat1 = sd_prop_diag1*corr_prop_mat1*sd_prop_diag1;

sd_prop_diag2 = diag([sd_prop_mean2,sd_prop_sd2,sd_prop_beta2]);
theta_prop_cov_mat2 = sd_prop_diag2*corr_prop_mat2*sd_prop_diag2;

sd_prop_diag3 = diag([sd_prop_mean3,sd_prop_sd3,sd_prop_beta3]);
theta_prop_cov_mat3 = sd_prop_diag3*corr_prop_mat3*sd_prop_diag3;

% Initial values of model parameters

theta_init = [(5),(5),(2)];

% Standard deviation of proposal distributions for infection times of
% symptomatic and asymptomatic infected hosts

t_i_prop_sd_symp1 = 8;
t_i_prop_sd_asymp1 = 15;
t_i_prop_sd_symp2 = 8;
t_i_prop_sd_asymp2 = 15;
t_i_prop_sd_symp3 = 7;
t_i_prop_sd_asymp3 = 8;

% Function handle giving the prior density for parameters theta

prior_fun = @(theta)prior_fun_indep(theta);

% Function handles used to update either the model parameters or augmented
% data during the fitting procedure

update_theta1 = @(theta,data_struct_augmented,ll_household)update_theta_fun(theta,data_struct_augmented,ll_household,ll_household_form,theta_prop_cov_mat1,prior_fun);
update_theta2 = @(theta,data_struct_augmented,ll_household)update_theta_fun(theta,data_struct_augmented,ll_household,ll_household_form,theta_prop_cov_mat2,prior_fun);
update_theta3 = @(theta,data_struct_augmented,ll_household)update_theta_fun(theta,data_struct_augmented,ll_household,ll_household_form,theta_prop_cov_mat3,prior_fun);
update_infection1 = @(theta,data_struct_augmented,ll_household)update_infection_fun_indep(theta,data_struct_augmented,ll_household,ll_household_form,t_i_prop_sd_symp1);
update_infection2 = @(theta,data_struct_augmented,ll_household)update_infection_fun_indep(theta,data_struct_augmented,ll_household,ll_household_form,t_i_prop_sd_symp2);
update_infection3 = @(theta,data_struct_augmented,ll_household)update_infection_fun_indep(theta,data_struct_augmented,ll_household,ll_household_form,t_i_prop_sd_symp3);
update_onset = @(theta,data_struct_augmented,ll_household)update_onset_fun(theta,data_struct_augmented,ll_household,ll_household_form);
update_asymp1 = @(theta,data_struct_augmented,ll_household)update_asymp_fun_indep(theta,data_struct_augmented,ll_household,ll_household_form,t_i_prop_sd_asymp1);
update_asymp2 = @(theta,data_struct_augmented,ll_household)update_asymp_fun_indep(theta,data_struct_augmented,ll_household,ll_household_form,t_i_prop_sd_asymp2);
update_asymp3 = @(theta,data_struct_augmented,ll_household)update_asymp_fun_indep(theta,data_struct_augmented,ll_household,ll_household_form,t_i_prop_sd_asymp3);

% Using data from first sub-interval, initialise augmented data and
% calculate the initial likelihood

data_struct_augmented_init1 = initialise_augmented_data_indep(data_struct_observed1);
ll_household_init1 = ll_household_form(theta_init,data_struct_augmented_init1);

% Run main parameter fitting procedure

tic
plotting = false; %set to true to plot output of fitting procedure as it is carried out (updating the plot once for each 1% of steps completed)
[theta_mat1,ll_vec1,output1] = fit_params(no_steps,steps_keep,update_theta1,update_infection1,update_onset,update_asymp1,empirical_summary_form,theta_init,data_struct_augmented_init1,ll_household_init1,plotting);
toc

% Repeat for data from other two sub-intervals

data_struct_augmented_init2 = initialise_augmented_data_indep(data_struct_observed2);
ll_household_init2 = ll_household_form(theta_init,data_struct_augmented_init2);

tic
[theta_mat2,ll_vec2,output2] = fit_params(no_steps,steps_keep,update_theta2,update_infection2,update_onset,update_asymp2,empirical_summary_form,theta_init,data_struct_augmented_init2,ll_household_init2,plotting);
toc

data_struct_augmented_init3 = initialise_augmented_data_indep(data_struct_observed3);
ll_household_init3 = ll_household_form(theta_init,data_struct_augmented_init3);

tic
[theta_mat3,ll_vec3,output3] = fit_params(no_steps,steps_keep,update_theta3,update_infection3,update_onset,update_asymp3,empirical_summary_form,theta_init,data_struct_augmented_init3,ll_household_init3,plotting);
toc

% Save results

save('../../Results/Split_data/param_fit_indep_splitdata.mat','theta_mat1','ll_vec1','output1','theta_mat2','ll_vec2','output2','theta_mat3','ll_vec3','output3')

rmpath('../../Data')
rmpath('../../Functions/MCMC')
rmpath('../../Functions/Indep')