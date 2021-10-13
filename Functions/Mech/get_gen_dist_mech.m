function f_gen_cheb = get_gen_dist_mech(params,t_range)

    % Calculate the generation time distribution for our mechanistic
    % approach with parameters given by params.
    
    % This function requires the Chebfun package to execute (freely
    % available at https://www.chebfun.org/download/).

    gamma = params(1); mu = params(2);
    k_inc = params(3); k_E = params(4); k_I = params(5);
    alpha = params(6);
    k_P = k_inc-k_E;

    C = k_inc*gamma*mu/(alpha*k_P*mu+k_inc*gamma);

    f_E = @(t)gampdf(t,k_E,1/(k_inc*gamma));
    f_E_cheb = chebfun(f_E,t_range);
    f_P = @(t)gampdf(t,k_P,1/(k_inc*gamma));
    f_P_cheb = chebfun(f_P,t_range);
    F_P = @(t)gamcdf(t,k_P,1/(k_inc*gamma));
    F_P_cheb = chebfun(F_P,t_range);
    F_I = @(t)gamcdf(t,k_I,1/(k_I*mu));
    F_I_cheb = chebfun(F_I,t_range);

    f1 = conv(f_P_cheb,F_I_cheb,'same','old');
    f_star_cheb1 = alpha*C*(1-F_P_cheb)+C*(F_P_cheb-f1);
    f_star_cheb = chebfun(@(t)f_star_cheb1(t).*(t>0),t_range);

    f_gen_cheb = conv(f_E_cheb,f_star_cheb,'same','old');

end