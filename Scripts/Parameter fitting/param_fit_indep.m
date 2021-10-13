% Fit the parameters of the independent transmission and symptoms model
% using data augmentation MCMC

% The vector of unknown model parameters, theta  = [m_gen,s_gen,beta_0], is
% estimated in the model fitting procedure. Here, m_gen and s_gen are the
% mean and standard deviation of the generation time distribution,
% respectively, and beta_0 is the overall infectiousness parameter

clear all; close all; clc;

addpath('../../Data')
addpath('../../Functions/MCMC')
addpath('../../Functions/Indep')

% Initialise random number generator

rng(3)

% Load data

load('../../Data/data.mat','data_struct_observed')

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

sd_prop_mean = 0.3;
sd_prop_sd = 2.5*sd_prop_mean;
sd_prop_beta = 0.25*sd_prop_mean;

% Correlation between proposal distributions

corr_prop_mat = eye(3);

% Covariance matrix of multivariate normal proposal distribution for the
% vector of fitted model parameters, theta

sd_prop_diag = diag([sd_prop_mean,sd_prop_sd,sd_prop_beta]);
theta_prop_cov_mat = sd_prop_diag*corr_prop_mat*sd_prop_diag;

% Initial values of model parameters

theta_init = [5,5,2];

% Standard deviation of proposal distributions for infection times of
% symptomatic and asymptomatic infected hosts

t_i_prop_sd_symp = 8;
t_i_prop_sd_asymp = 13;

% Function handle giving the prior density for parameters theta

prior_fun = @(theta)prior_fun_indep(theta);

% Function handles used to update either the model parameters or augmented
% data during the fitting procedure

update_theta = @(theta,data_struct_augmented,ll_household)update_theta_fun(theta,data_struct_augmented,ll_household,ll_household_form,theta_prop_cov_mat,prior_fun);
update_infection = @(theta,data_struct_augmented,ll_household)update_infection_fun_indep(theta,data_struct_augmented,ll_household,ll_household_form,t_i_prop_sd_symp);
update_onset = @(theta,data_struct_augmented,ll_household)update_onset_fun(theta,data_struct_augmented,ll_household,ll_household_form);
update_asymp = @(theta,data_struct_augmented,ll_household)update_asymp_fun_indep(theta,data_struct_augmented,ll_household,ll_household_form,t_i_prop_sd_asymp);

% Initialise augmented data and calculate the initial likelihood

data_struct_augmented_init = initialise_augmented_data_indep(data_struct_observed);
ll_household_init = ll_household_form(theta_init,data_struct_augmented_init);

% Run main parameter fitting procedure

tic
plotting = false; %set to true to plot output of fitting procedure as it is carried out (updating the plot once for each 1% of steps completed)
[theta_mat,ll_vec,output] = fit_params(no_steps,steps_keep,update_theta,update_infection,update_onset,update_asymp,empirical_summary_form,theta_init,data_struct_augmented_init,ll_household_init,plotting);
toc

% Save results

save('../../Results/param_fit_indep.mat','theta_mat','ll_vec','output')

rmpath('../../Data')
rmpath('../../Functions/MCMC')
rmpath('../../Functions/Indep')