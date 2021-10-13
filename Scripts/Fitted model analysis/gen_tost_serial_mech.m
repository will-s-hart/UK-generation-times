% Calculate distributions of the generation time, time from onset of
% symptoms to transmission (TOST) and serial interval for the mechanistic
% model, using point estimate (posterior mean) parameter values

% This script requires the Chebfun package to run (freely available at
% https://www.chebfun.org/download/)

clear all; close all; clc;
splitting on

addpath('../../Data')
addpath('../../Results')
addpath('../../Functions/Mech')

% Load point estimate parameter values

load('../../Results/mcmc_posterior_mech.mat','params_best')

% Generation time

t_range = [-100,-50,-25,0,25,50,100];
f_gen_mech = get_gen_dist_mech(params_best,t_range);

% TOST

f_tost_mech = chebfun(@(t)f_tost_form_mech(t,params_best),t_range);

% Serial interval

f_serial_mech = get_serial_dist_mech(params_best,t_range);

% Save results

save('../../Results/gen_tost_serial_mech.mat','f_gen_mech','f_tost_mech','f_serial_mech')

rmpath('../../Data')
rmpath('../../Results')
rmpath('../../Functions/Mech')