% Import household transmission data from the file
% "Supplementary_Data.xlsx", and apply preliminary formatting steps.

clear all; close all; clc;

% Import data

T = readtable('Supplementary_Data.xlsx');

% Extract the household number, household size and symptom onset date for
% each host (day zero taken to be date on which index case swabbed)

household_no_all = T.hoconumber;
household_size_all = T.householdsize;
d_s_all = T.case_swab_ill;

month_all = T.month_case_swab;

% Logical vectors indicating whether each host
% infected/uninfected/inconclusive, asymptomatic/symptomatic

inf_all = (T.infected==1);
uninf_all = (T.infected==0);
inconclusive_all = ~(inf_all|uninf_all);
symp_all = (inf_all&(T.symptoms==1));
asymp_all = (inf_all&(T.symptoms==0));

% Uninfected or asymptomatic individuals considered to develop symptoms at
% time infinity

d_s_all(asymp_all) = inf;
d_s_all(uninf_all) = inf;

% Find index cases

status_all = T.status;
l = cellfun(@length,status_all);
index_all = (l==4);

% Obtain right bounds for infection date of each host

d_iR_all = inf*ones(length(d_s_all),1);
d_iR_all(symp_all) = d_s_all(symp_all); %can't be infected after day of onset
d_iR_all(index_all) = min(d_iR_all(index_all),0); %since date where index swabbed positive taken as time zero

swab1_positive = (T.swab1==1);
swab2_positive = (T.swab2==1);

d_iR_all(swab1_positive) = min(d_iR_all(swab1_positive),T.case_swab_swab1(swab1_positive)); %must have been infected before positive test
d_iR_all(swab2_positive) = min(d_iR_all(swab2_positive),T.case_swab_swab2(swab2_positive)); %must have been infected before positive test

% Discard data from hosts not known to have been infected or not

keep = symp_all|asymp_all|uninf_all;

household_no_all = household_no_all(keep);
household_size_all = household_size_all(keep); %includes thrown out hosts
month_all = month_all(keep);
d_s_all = d_s_all(keep);
d_iR_all = d_iR_all(keep);
symp_all = symp_all(keep);
asymp_all = asymp_all(keep);
uninf_all = uninf_all(keep);
inf_all = ~uninf_all;
index_all = index_all(keep);

% Order data in each household by symptom onset date, also putting
% asymptomatic hosts before uninfected hosts

[~,order] = sortrows([household_no_all,uninf_all,d_s_all]);

household_no_all = household_no_all(order);
household_size_all = household_size_all(order);
month_all = month_all(order);
d_s_all = d_s_all(order);
d_iR_all = d_iR_all(order);
symp_all = symp_all(order);
asymp_all = asymp_all(order);
uninf_all = uninf_all(order);
inf_all = inf_all(order);
index_all = index_all(order);

% Households kept once data from asymptomatic/inconclusive hosts thrown out

households_vec_all = unique(household_no_all);
no_households_all = length(households_vec_all);

% Assumed maximum possible gap between successive symptom onset dates
% within a household

max_onsetdiff = 28;

% Separate out clusters when the maximum onset difference is exceeded

no_households_incl = 0;
households_incl_old = [];
households_incl_new = [];
household_no_new_incl = [];
household_sizes_old_incl = []; %counting thrown out hosts of unknown status
household_sizes_new_incl = []; %not counting thrown out hosts
household_size_indiv_old_incl = [];
household_size_indiv_new_incl = [];
household_months_incl = [];
d_s_incl = [];
d_iR_incl = [];
symp_incl = [];
asymp_incl = [];
uninf_incl = [];
onset_diffs_incl = [];

for household = households_vec_all'
    
    % Extract indices, symptom onset dates and status of all individuals in
    % household
    
    household_inds = find(household_no_all==household);
    household_d_s = d_s_all(household_inds);
    household_d_iR = d_iR_all(household_inds);
    household_symp = symp_all(household_inds);
    household_asymp = asymp_all(household_inds);
    household_uninf = uninf_all(household_inds);
        
    % Household size, number who developed symptoms, month in which index
    % case was recruited
    
    no_in_household_old = household_size_all(household_inds(1)); %counting thrown out hosts of unknown status
    no_in_household = length(household_inds); %not counting thrown out hosts
    no_symp_in_household = sum(household_symp);
    month_household = month_all(household_inds(1));
    
    % Only keep the data from the household if (i) it contains more than
    % one household member of known infection status, and (ii) each gap
    % between successive symptom onset dates is no larger than the
    % permitted maximum
    
    household_onset_diffs = diff(household_d_s(household_symp));
    household_onset_diff_max = max(household_onset_diffs);
    
    if (no_in_household>1)&&((isempty(household_onset_diff_max)||household_onset_diff_max<=max_onsetdiff))
        no_households_incl = no_households_incl + 1;
        households_incl_old = [households_incl_old;household];
        households_incl_new = [households_incl_new;no_households_incl];
        household_no_new_incl = [household_no_new_incl;no_households_incl*ones(no_in_household,1)];
        household_sizes_old_incl = [household_sizes_old_incl;no_in_household_old];
        household_sizes_new_incl = [household_sizes_new_incl;no_in_household];
        household_size_indiv_old_incl = [household_size_indiv_old_incl;no_in_household_old*ones(no_in_household,1)];
        household_size_indiv_new_incl = [household_size_indiv_new_incl;no_in_household*ones(no_in_household,1)];
        household_months_incl = [household_months_incl;month_household];
        d_s_incl = [d_s_incl;household_d_s];
        d_iR_incl = [d_iR_incl;household_d_iR];
        symp_incl = [symp_incl;household_symp];
        asymp_incl = [asymp_incl;household_asymp];
        uninf_incl = [uninf_incl;household_uninf];
        onset_diffs_incl = [onset_diffs_incl;household_onset_diffs];
    end
end

% Convert vectors indicating infection/symptom status to logical vectors

symp_incl = logical(symp_incl);
asymp_incl = logical(asymp_incl);
uninf_incl = logical(uninf_incl);

% Save data

save('data_initial.mat','household_no_new_incl','household_sizes_old_incl','household_sizes_new_incl','household_size_indiv_old_incl','household_size_indiv_new_incl','household_months_incl','d_s_incl','d_iR_incl','symp_incl','asymp_incl','uninf_incl')