% Fit the parameters of the mechanistic model to the data from three
% sub-intervals of the study interval (the suffixes 1,2,3 in variable names
% refer to the different time-intervals), using data augmentation MCMC

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

load('../../Data/data_split.mat','data_split')
data_struct_observed1 = data_split{1};
data_struct_observed2 = data_split{2};
data_struct_observed3 = data_split{3};

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

sd_prop_p_E1 = 0.085;
sd_prop_mu_inv1 = 10*sd_prop_p_E1;
sd_prop_alpha1 = 10*sd_prop_p_E1;
sd_prop_beta1 = 1.4*sd_prop_p_E1;

% Correlation between proposal distributions

corr_p_E_mu_inv1 = 0;
corr_p_E_alpha1 = 0.5;
corr_p_E_beta1 = 0;
corr_mu_inv_alpha1 = 0.5;
corr_mu_inv_beta1 = 0;
corr_alpha_beta1 = 0;

corr_prop_mat1 = [1,corr_p_E_mu_inv1,corr_p_E_alpha1,corr_p_E_beta1;corr_p_E_mu_inv1,1,corr_mu_inv_alpha1,corr_mu_inv_beta1;corr_p_E_alpha1,corr_mu_inv_alpha1,1,corr_alpha_beta1;corr_p_E_beta1,corr_mu_inv_beta1,corr_alpha_beta1,1];

% Covariance matrix of multivariate normal proposal distribution for the
% vector of fitted model parameters, theta

sd_prop_diag1 = diag([sd_prop_p_E1,sd_prop_mu_inv1,sd_prop_alpha1,sd_prop_beta1]);
theta_prop_cov_mat1 = sd_prop_diag1*corr_prop_mat1*sd_prop_diag1;

theta_prop_cov_mat2 = (1.3^2)*theta_prop_cov_mat1;
theta_prop_cov_mat3 = (1.25^2)*theta_prop_cov_mat1;

% Initial values of model parameters

theta_init = [0.5,1/(0.18),(3.5),(2)];

% Standard deviation of proposal distributions for infection times of
% infected hosts, and for simultaneous updates of the times at which
% asymptomatic hosts become infected and enter the I stage (holding the
% difference between the two times constant; note the distinction between
% the P and I stages has no epidemiological meaning for asymptomatic hosts,
% but the observed data are augmented with the time of entering the I stage
% when fitting the model for convenience)

t_i_prop_sd1 = 9;
t_prop_sd_asymp1 = 20;
t_i_prop_sd2 = 9;
t_prop_sd_asymp2 = 20;
t_i_prop_sd3 = 8;
t_prop_sd_asymp3 = 16;

% Function handle giving the prior density for parameters theta

prior_fun = @(theta)prior_fun_mech(theta);

% Function handles used to update either the model parameters or augmented
% data during the fitting procedure

update_theta1 = @(theta,data_struct_augmented,ll_household)update_theta_fun(theta,data_struct_augmented,ll_household,ll_household_form,theta_prop_cov_mat1,prior_fun);
update_theta2 = @(theta,data_struct_augmented,ll_household)update_theta_fun(theta,data_struct_augmented,ll_household,ll_household_form,theta_prop_cov_mat2,prior_fun);
update_theta3 = @(theta,data_struct_augmented,ll_household)update_theta_fun(theta,data_struct_augmented,ll_household,ll_household_form,theta_prop_cov_mat3,prior_fun);
update_infection1 = @(theta,data_struct_augmented,ll_household)update_infection_fun_mech(theta,data_struct_augmented,ll_household,ll_household_form,t_i_prop_sd1);
update_infection2 = @(theta,data_struct_augmented,ll_household)update_infection_fun_mech(theta,data_struct_augmented,ll_household,ll_household_form,t_i_prop_sd2);
update_infection3 = @(theta,data_struct_augmented,ll_household)update_infection_fun_mech(theta,data_struct_augmented,ll_household,ll_household_form,t_i_prop_sd3);
update_onset = @(theta,data_struct_augmented,ll_household)update_onset_fun(theta,data_struct_augmented,ll_household,ll_household_form);
update_asymp1 = @(theta,data_struct_augmented,ll_household)update_asymp_fun_mech(theta,data_struct_augmented,ll_household,ll_household_form,t_prop_sd_asymp1);
update_asymp2 = @(theta,data_struct_augmented,ll_household)update_asymp_fun_mech(theta,data_struct_augmented,ll_household,ll_household_form,t_prop_sd_asymp2);
update_asymp3 = @(theta,data_struct_augmented,ll_household)update_asymp_fun_mech(theta,data_struct_augmented,ll_household,ll_household_form,t_prop_sd_asymp3);

% Using data from first sub-interval, initialise augmented data and
% calculate the initial likelihood

data_struct_augmented_init1 = initialise_augmented_data_mech(data_struct_observed1);
ll_household_init1 = ll_household_form(theta_init,data_struct_augmented_init1);

% Run main parameter fitting procedure

tic
plotting = false; %set to true to plot output of fitting procedure as it is carried out (updating the plot once for each 1% of steps completed)
[theta_mat1,ll_vec1,output1] = fit_params(no_steps,steps_keep,update_theta1,update_infection1,update_onset,update_asymp1,empirical_summary_form,theta_init,data_struct_augmented_init1,ll_household_init1,plotting);
toc

% Repeat for data from other two sub-intervals

data_struct_augmented_init2 = initialise_augmented_data_mech(data_struct_observed2);
ll_household_init2 = ll_household_form(theta_init,data_struct_augmented_init2);

tic
[theta_mat2,ll_vec2,output2] = fit_params(no_steps,steps_keep,update_theta2,update_infection2,update_onset,update_asymp2,empirical_summary_form,theta_init,data_struct_augmented_init2,ll_household_init2,plotting);
toc

data_struct_augmented_init3 = initialise_augmented_data_mech(data_struct_observed3);
ll_household_init3 = ll_household_form(theta_init,data_struct_augmented_init3);

tic
[theta_mat3,ll_vec3,output3] = fit_params(no_steps,steps_keep,update_theta3,update_infection3,update_onset,update_asymp3,empirical_summary_form,theta_init,data_struct_augmented_init3,ll_household_init3,plotting);
toc

% Save results

save('../../Results/Split_data/param_fit_mech_splitdata.mat','theta_mat1','ll_vec1','output1','theta_mat2','ll_vec2','output2','theta_mat3','ll_vec3','output3')

rmpath('../../Data')
rmpath('../../Functions/MCMC')
rmpath('../../Functions/Mech')