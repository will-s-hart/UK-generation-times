function f_serial_cheb = get_serial_dist_mech(params,t_range)

    % Calculate the generation time distribution for our mechanistic
    % approach with parameters given by params.
    
    % This function requires the Chebfun package to execute (freely
    % available at https://www.chebfun.org/download/).

    gamma = params(1);
    k_inc = params(3);

    f_inc = @(t_inc)gampdf(t_inc,k_inc,1/(k_inc*gamma));
    f_inc_cheb = chebfun(f_inc,t_range);

    f_tost = @(t)f_tost_form_mech(t,params);
    f_tost_cheb = chebfun(f_tost,t_range);

    f_serial_cheb = conv(f_inc_cheb,f_tost_cheb,'same','old');

end