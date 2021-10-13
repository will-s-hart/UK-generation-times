% Separately format data from three time intervals within the study period
% into structure arrays.

% Throughout, arrays with the suffix "_dir" are ordered according to (i)
% the household number assigned to each household, and (ii) the (unknown)
% order in which household members became infected. Otherwise, arrays are
% ordered according to (i) household number, and (ii) the (known) order in
% which household members developed symptoms.

clear all; close all; clc;

% Load imported data in file "data_initial.mat" (this file is created by
% running "import_data.m")

load('data_initial.mat','household_no_new_incl','household_sizes_old_incl','household_sizes_new_incl','household_size_indiv_old_incl','household_size_indiv_new_incl','household_months_incl','d_s_incl','d_iR_incl','uninf_incl','symp_incl','asymp_incl')

% Month in which each household was recruited

household_months = household_months_incl;

% Divide households and hosts into three subintervals

cutoff1 = 4.5;
cutoff2 = 8.5;
households1 = find(household_months<cutoff1);
households2 = find((household_months>=cutoff1)&(household_months<cutoff2));
households3 = find(household_months>=cutoff2);
hosts1 = find(household_months(household_no_new_incl)<cutoff1);
hosts2 = find((household_months(household_no_new_incl)>=cutoff1)&(household_months(household_no_new_incl)<cutoff2));
hosts3 = find(household_months(household_no_new_incl)>=cutoff2);

% Cell array with households and hosts in each subinterval

household_split = {households1,households2,households3};
host_split = {hosts1,hosts2,hosts3};

% Cell array to populate with data from each subinterval

data_split = cell(1,3);

for k = 1:3
    
    % Hosts and households within subinterval
    
    households_incl = household_split{k};
    hosts_incl = host_split{k};
    
    % Number of hosts and households
    
    no_hosts = length(hosts_incl);
    no_households = length(households_incl);
    
    % Extract data from hosts and households within subinterval
    
    household_sizes_incl = household_sizes_new_incl(households_incl);
    household_sizes_full = household_sizes_old_incl(households_incl);
    household_size_indiv_incl = household_size_indiv_new_incl(hosts_incl);
    household_size_indiv_full = household_size_indiv_old_incl(hosts_incl);
    household_no = household_no_new_incl(hosts_incl);
    
    d_s = d_s_incl(hosts_incl);
    d_iR = d_iR_incl(hosts_incl);
    uninfected_dir = uninf_incl(hosts_incl);
    infected_dir = ~uninfected_dir;
    symp = symp_incl(hosts_incl);
    asymp = asymp_incl(hosts_incl);
    
    % Relabel households within subinterval
    
    households_incl_new = (1:no_households)';
    
    household_no_new = zeros(no_households,1);
    for i = 1:no_hosts
        
        household_no_new(i) = households_incl_new(households_incl==household_no(i));
    end
    household_no = household_no_new;
    
    % Lower and upper bounds for symptom onset and infection times

    t_sL = d_s-0.5;
    t_sR = d_s+0.5;

    t_iL = -inf*ones(no_hosts,1);
    t_iR = d_iR+0.5;
    
    % Indicator matrix showing which household each individual belongs to
    
    blocks = mat2cell(sparse(ones(no_hosts,1)),household_sizes_incl,1);
    household_indicator_mat = blkdiag(blocks{:});
    
    % Calculate the number of infected and symptomatic/asymptomatic hosts
    % in each household, and working in the (unknown) order of infection,
    % determine the possible (household) infectors for each individual
    % (i.e., household members who were infected before that individual).
    
    no_infected_in_household = zeros(no_households,1);
    no_symp_in_household = zeros(no_households,1);
    no_asymp_in_household = zeros(no_households,1);
    no_poss_infectors_dir = zeros(no_hosts,1);
    poss_infectors_dir_cell = cell(no_hosts,1);
    primary_dir = false(no_hosts,1);

    for i = 1:no_households

        in_household = (household_no==i);
        infected_hosts_in_household = find(in_household.*infected_dir);
        symp_hosts_in_household = find(in_household.*symp);
        asymp_hosts_in_household = find(in_household.*asymp);
        uninfected_hosts_in_household = find(in_household.*uninfected_dir);

        no_infected_in_household(i) = length(infected_hosts_in_household);
        no_symp_in_household(i) = length(symp_hosts_in_household);
        no_asymp_in_household(i) = length(asymp_hosts_in_household);

        no_poss_infectors_dir(infected_hosts_in_household(1)) = 1;
        poss_infectors_dir_cell{infected_hosts_in_household(1)} = 0; %use 0 to denote infection from outside the household
        primary_dir(infected_hosts_in_household(1)) = true;

        for j = 2:length(infected_hosts_in_household)
            poss_infectors_dir_cell{infected_hosts_in_household(j)} = infected_hosts_in_household(1:(j-1));
            no_poss_infectors_dir(infected_hosts_in_household(j)) = length(poss_infectors_dir_cell{infected_hosts_in_household(j)});
        end

        for j = 1:length(uninfected_hosts_in_household)
            poss_infectors_dir_cell{uninfected_hosts_in_household(j)} = infected_hosts_in_household;
            no_poss_infectors_dir(uninfected_hosts_in_household(j)) = length(infected_hosts_in_household);
        end
    end
    
    % Logical arrays indicating whether or not each household contained any
    % symptomatic or asymptomatic infected hosts
    
    symp_in_household = any(no_symp_in_household,2);
    asymp_in_household = any(no_asymp_in_household,2);    
    
    % Create a vector v, enumerating each values within the cell array
    % listing the possible infectors for each individual (in the order of
    % infection)
    
    v = cell2mat(poss_infectors_dir_cell);
    
    % Indicator matrix giving the index of the potential infector
    % corresponding to each entry of v
    
    M1 = spalloc(length(v),no_hosts,length(v)-no_households);

    for i = 1:length(v)

        j = v(i);

        if j > 0
            M1(i,j) = 1;
        end
    end

    % Indicator matrix giving the index of the potential infectee
    % corresponding to each entry of v

    blocks = mat2cell(sparse(ones(length(v),1)),no_poss_infectors_dir,1);
    M2 = blkdiag(blocks{:}); %indicator of recipient corresponding to each entry of v
    
    % Logical vectors indicating, respectively, whether or not each entry of v
    % corresponds to: (i) the potential infector of a host who actually became
    % infected; (ii) the infection of the primary case (denoted by a 0 in v);
    % (iii) the potential infector of a non-primary case who became infected
    
    v_to_infected_indicator = logical(M2*infected_dir);
    v_to_primary_indicator = (v==0);
    v_to_recipient_indicator = (v_to_infected_indicator&~v_to_primary_indicator);

    % The size of the household that each entry of v corresponds to
    
    household_size_v = M2*household_size_indiv_full;
    
    % Structure array containing information about possible infectors
    
    poss_infectors_dir.cell = poss_infectors_dir_cell;
    poss_infectors_dir.all = v;
    poss_infectors_dir.from_indicator_mat = M1;
    poss_infectors_dir.to_indicator_mat = M2;
    poss_infectors_dir.to_primary_indicator = v_to_primary_indicator;
    poss_infectors_dir.to_infected_indicator = v_to_infected_indicator;
    poss_infectors_dir.to_recipient_indicator = v_to_recipient_indicator;
    poss_infectors_dir.household_size = household_size_v;

    % Structure array containing all observed data from the subinterval

    data_struct_observed.t_iL = t_iL;
    data_struct_observed.t_iR = t_iR;
    data_struct_observed.t_sL = t_sL;
    data_struct_observed.t_sR = t_sR;

    data_struct_observed.household_sizes_incl = household_sizes_incl;
    data_struct_observed.household_sizes_full = household_sizes_full;
    data_struct_observed.household_size_indiv_incl = household_size_indiv_incl;
    data_struct_observed.household_size_indiv_full = household_size_indiv_full;
    data_struct_observed.household_no = household_no;
    data_struct_observed.household_indicator_mat = household_indicator_mat;
    data_struct_observed.no_infected_in_household = no_infected_in_household;
    data_struct_observed.no_symp_in_household = no_symp_in_household;
    data_struct_observed.symp_in_household = symp_in_household;
    data_struct_observed.no_asymp_in_household = no_asymp_in_household;
    data_struct_observed.asymp_in_household = asymp_in_household;
    data_struct_observed.primary_dir = primary_dir;
    data_struct_observed.infected_dir = infected_dir;
    data_struct_observed.symp = symp;
    data_struct_observed.asymp = asymp;
    data_struct_observed.no_poss_infectors_dir = no_poss_infectors_dir;
    data_struct_observed.poss_infectors_dir = poss_infectors_dir;   
    
    % Populate cell array with data from the subinterval
    
    data_split{k} = data_struct_observed;
end

% Save data to .mat file

save('data_split.mat','data_split')