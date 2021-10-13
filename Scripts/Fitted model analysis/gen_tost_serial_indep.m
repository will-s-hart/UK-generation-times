% Calculate distributions of the generation time, time from onset of
% symptoms to transmission (TOST) and serial interval for the independent
% transmission and symptoms model, using point estimate (posterior mean)
% parameter values

% This script requires the Chebfun package to run (freely available at
% https://www.chebfun.org/download/)

clear all; close all; clc;
splitting on

addpath('../../Data')
addpath('../../Results')

% Load incubation period distribution

load('../../Data/assumed_parameters.mat','f_inc_logn')
f_inc = f_inc_logn;

% Load point estimate parameter values

load('../../Results/mcmc_posterior_indep.mat','mean_best','sd_best')

logn_mu = @(m,s)log(m^2/sqrt(s^2+m^2));
logn_sigma = @(m,s)sqrt(log(1+s^2/m^2));

mu_best = logn_mu(mean_best,sd_best);
sigma_best = logn_sigma(mean_best,sd_best);

% Generation time

t_range = [-100,-50,-25,0,25,50,100];
f_gen_indep = @(t,params)lognpdf(t,mu_best,sigma_best);
f_gen_indep = chebfun(f_gen_indep,t_range);

% TOST

f_inc = chebfun(f_inc,t_range);
f_tost_indep = conv(f_gen_indep,flipud(f_inc),'same','old');

% Serial interval

f_tost_indep = merge(f_tost_indep);
f_serial_indep = conv(f_inc,f_tost_indep,'same','old');

% Save results

save('../../Results/gen_tost_serial_indep.mat','f_gen_indep','f_tost_indep','f_serial_indep')

rmpath('../../Data')
rmpath('../../Results')