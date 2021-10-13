function empirical_summary = empirical_summary_indep(beta0,rho,x_A,f_gen,data_struct_augmented)

    % Estimate the mean and standard deviation of realised household
    % generation times, in addition to the proportion of presymptomatic
    % transmissions, for the independent transmission and symptoms model
    % with augmented data data_struct_augmented.
    
    % Throughout, arrays with the suffix "_dir" are ordered according to
    % (i) the household number assigned to each household, and (ii) the
    % (unknown) order in which household members became infected.

    % Augmented infection and symptom onset times

    t_i_dir = data_struct_augmented.t_i_dir;
    t_s_dir = data_struct_augmented.t_s_dir;
    
    % Logical vectors indicating symptom status
    
    symp_dir = data_struct_augmented.symp_dir;
    asymp_dir = data_struct_augmented.asymp_dir;
    
    % Information about possible infectors
    
    poss_infectors_dir = data_struct_augmented.poss_infectors_dir;
    v = poss_infectors_dir.all; %list of all possible infectors for each host in order
    M_from = poss_infectors_dir.from_indicator_mat;
    M_to = poss_infectors_dir.to_indicator_mat;
    to_primary_indicator_v = poss_infectors_dir.to_primary_indicator;
    to_recipient_indicator_v = poss_infectors_dir.to_recipient_indicator;
    household_size_v = poss_infectors_dir.household_size;

    from_symp_indicator_v = logical(M_from*symp_dir);
    from_asymp_indicator_v = logical(M_from*asymp_dir);

    beta_v = beta0./(household_size_v.^rho);
    beta_v(to_primary_indicator_v) = 0;
    beta_v(from_asymp_indicator_v) = x_A*beta_v(from_asymp_indicator_v);
    
    % Vectors of every possible generation time for each infectee (with
    % entries corresponding to infection by every possible infector), and
    % of the incubation period (which is infinite for asymptomatic hosts)
    % of the corresponding infectee.
    
    t_gen_contribs = (M_to-M_from)*t_i_dir;
    t_inc_contribs = M_from*(t_s_dir-t_i_dir);
    
    % Calculate likelihood contributions corresponding to infection of each
    % individual
    
    t_gen_recipient_contribs = t_gen_contribs(to_recipient_indicator_v);
    beta_recipient_contribs = beta_v(to_recipient_indicator_v);

    L2a_contribs = zeros(length(v),1);
    L2a_contribs(to_recipient_indicator_v) = beta_recipient_contribs.*f_gen(t_gen_recipient_contribs);
    L2a = M_to'*L2a_contribs;
    
    % Non-trivial possible generation times and times from symptom onset to
    % transmission (TOST)
    
    t_gen = t_gen_recipient_contribs;
    t_tost = t_gen_contribs(from_symp_indicator_v&to_recipient_indicator_v)-t_inc_contribs(from_symp_indicator_v&to_recipient_indicator_v);
    
    % Use likelihood contributions to calculate weights indicating the
    % relative probabilities of infection by different possible infectiors
    
    weights_all = L2a_contribs./(M_to*L2a);
    weights_gen = weights_all(to_recipient_indicator_v);
    weights_tost = weights_all(from_symp_indicator_v&to_recipient_indicator_v);

    weights_gen = weights_gen/sum(weights_gen);
    weights_tost = weights_tost/sum(weights_tost);
    
    % Calculate the mean and standard deviation of realised generation
    % times, using the computed weights to account for different
    % probabilities that individuals were infected by different infectors.
    
    m = sum(t_gen.*weights_gen);
    s = sqrt(sum(t_gen.*t_gen.*weights_gen)-m^2);
    p = sum((t_tost<0).*weights_tost);
    
    empirical_summary = [m,s,p];
end