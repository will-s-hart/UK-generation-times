% Fit the parameters of the mechanistic model using data augmentation MCMC

% The vector of unknown model parameters, theta  =
% [p_E,mu_inv,alpha,beta_0], is estimated in the model fitting procedure.
% Here, p_E=k_E/k_inc is the ratio between the shape parameters of the
% Gamma distributed latent (E) and incubation (combined E and P) periods,
% mu_inv=1/mu is the mean symptomatic infectious (I) period, alpha=alpha_P
% is the ratio of the transmission rates during the presymptomatic
% infectious (P) and symptomatic infectious (I) periods, and beta_0 is the
% overall infectiousness parameter

clear all; close all; clc;

addpath('../../Data')
addpath('../../Functions/MCMC')
addpath('../../Functions/Mech')

% Initialise random number generator

rng(10)

% Load data

load('../../Data/data.mat','data_struct_observed')

% Load incubation period distribution function and assumed parameters

load('../../Data/assumed_parameters.mat','f_inc_gam','params_known','k_inc','gamma')
f_inc = f_inc_gam;

% Number of steps in chain, keeping only some steps after burn-in and
% thinning

no_steps = 10000000;
steps_keep = ((no_steps/5)+1):100:no_steps;

% Function handles giving (i) the expected infectiousness at each time
% since symptom onset x, conditional on incubation period t_inc; (ii) the
% integral of this infectiousness profile between times since onset
% (-infinity) and x; and (iii) the integral over all times since infection;
% for parameters theta

b_cond_form = @(x,t_inc,household_size,asymp,theta)b_cond_form_mech(x,t_inc,household_size,asymp,get_params_mech(theta,params_known));
B_cond_form = @(x,t_inc,household_size,asymp,theta)b_int_cond_form_mech(x,t_inc,household_size,asymp,get_params_mech(theta,params_known));
mean_transmissions_form = @(t_inc,household_size,asymp,theta)mean_transmissions_form_mech(t_inc,household_size,asymp,get_params_mech(theta,params_known));

% Function handle giving the log-likelihood for parameters theta and
% augmented data data_struct_augmented

ll_household_form = @(theta,data_struct_augmented)log_likelihood_household_mech(f_inc,@(x,t_inc,household_size,asymp)b_cond_form(x,t_inc,household_size,asymp,theta),@(x,t_inc,household_size,asymp)B_cond_form(x,t_inc,household_size,asymp,theta),@(t_inc,household_size,asymp)mean_transmissions_form(t_inc,household_size,asymp,theta),data_struct_augmented);

% Function handle giving the estimated mean and standard deviation of
% realised household generation times, and the presymptomatic proportion of
% transmissions, for parameters theta and augmented data
% data_struct_augmented

empirical_summary_form = @(theta,data_struct_augmented)empirical_summary_mech(@(x,t_inc,household_size,asymp)b_cond_form(x,t_inc,household_size,asymp,theta),data_struct_augmented);

% Standard deviations of proposal distributions for individual model
% parameters

sd_prop_p_E = 0.075;
sd_prop_mu_inv = 10*sd_prop_p_E;
sd_prop_alpha = 10*sd_prop_p_E;
sd_prop_beta = 1.4*sd_prop_p_E;

% Correlation between proposal distributions

corr_p_E_mu_inv = 0;
corr_p_E_alpha = 0.5;
corr_p_E_beta = 0;
corr_mu_inv_alpha = 0.5;
corr_mu_inv_beta = 0;
corr_alpha_beta = 0;

corr_prop_mat = [1,corr_p_E_mu_inv,corr_p_E_alpha,corr_p_E_beta;corr_p_E_mu_inv,1,corr_mu_inv_alpha,corr_mu_inv_beta;corr_p_E_alpha,corr_mu_inv_alpha,1,corr_alpha_beta;corr_p_E_beta,corr_mu_inv_beta,corr_alpha_beta,1];

% Covariance matrix of multivariate normal proposal distribution for the
% vector of fitted model parameters, theta

sd_prop_diag = diag([sd_prop_p_E,sd_prop_mu_inv,sd_prop_alpha,sd_prop_beta]);
theta_prop_cov_mat = sd_prop_diag*corr_prop_mat*sd_prop_diag;

% Initial values of model parameters

theta_init = [0.5,1/(0.18),(3.5),(2)];

% Standard deviation of proposal distributions for infection times of
% infected hosts, and for simultaneous updates of the times at which
% asymptomatic hosts become infected and enter the I stage (holding the
% difference between the two times constant; note the distinction between
% the P and I stages has no epidemiological meaning for asymptomatic hosts,
% but the observed data are augmented with the time of entering the I stage
% when fitting the model for convenience)

t_i_prop_sd = 9;
t_prop_sd_asymp = 18;

% Function handle giving the prior density for parameters theta

prior_fun = @(theta)prior_fun_mech(theta);

% Function handles used to update either the model parameters or augmented
% data during the fitting procedure

update_theta = @(theta,data_struct_augmented,ll_household)update_theta_fun(theta,data_struct_augmented,ll_household,ll_household_form,theta_prop_cov_mat,prior_fun);
update_infection = @(theta,data_struct_augmented,ll_household)update_infection_fun_mech(theta,data_struct_augmented,ll_household,ll_household_form,t_i_prop_sd);
update_onset = @(theta,data_struct_augmented,ll_household)update_onset_fun(theta,data_struct_augmented,ll_household,ll_household_form);
update_asymp = @(theta,data_struct_augmented,ll_household)update_asymp_fun_mech(theta,data_struct_augmented,ll_household,ll_household_form,t_prop_sd_asymp);

% Initialise augmented data and calculate the initial likelihood

data_struct_augmented_init = initialise_augmented_data_mech(data_struct_observed);
ll_household_init = ll_household_form(theta_init,data_struct_augmented_init);

% Run main parameter fitting procedure

tic
plotting = false; %set to true to plot output of fitting procedure as it is carried out (updating the plot once for each 1% of steps completed)
[theta_mat,ll_vec,output] = fit_params(no_steps,steps_keep,update_theta,update_infection,update_onset,update_asymp,empirical_summary_form,theta_init,data_struct_augmented_init,ll_household_init,plotting);
toc

% Save results

save('../../Results/param_fit_mech.mat','theta_mat','ll_vec','output')

rmpath('../../Data')
rmpath('../../Functions/MCMC')
rmpath('../../Functions/Mech')