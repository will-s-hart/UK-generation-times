function p_vec = get_presymp_trans_probs_indep_logn(gen_mu_vec,gen_sigma_vec,F_inc)
    
    % Numerically calculate the proportion of presymptomatic transmissions
    % for the independent transmission and symptoms model, for each set of
    % generation time parameters given by the entries in gen_mu_vec and
    % gen_sigma_vec, by integrating the generation time distribution
    % weighted by the proportion of hosts who have developed symptoms.
    
    % Grid of times since infection to integrate over
    
    t_max = 50;
    dt = 0.1;
    t_vec = ((dt/2):dt:(t_max-dt/2))';

    % Grids of pairs of values of the parameters of the generation time
    % distribution, and of the time since infection
    
    [gen_mu_grid,t_grid] = ndgrid(gen_mu_vec,t_vec);
    [gen_sigma_grid,~] = ndgrid(gen_sigma_vec,t_vec);

    % Evaluate the density of the generation time and the cumulative
    % distribution of the incubation period, at each grid point
    
    f_gen_grid = lognpdf(t_grid,gen_mu_grid,gen_sigma_grid);
    F_inc_vec = F_inc(t_vec);
    
    % Numerically integrate over the generation time distribution, weighted
    % by the proportion of hosts who have developed symptoms, to calculate
    % the proportion of presymptomatic transmissions for each parameter
    % set.
    
    p_vec = 1-f_gen_grid*F_inc_vec*dt;

end