function p = prior_fun_mech(theta)
    
    % Calculate the prior density for the independent transmission and
    % symptoms model with estimated parameters theta

    p_p_E = betapdf(theta(1),2.1,2.1); %median 0.5, 95% CI [0.1,0.9]
    p_mu_inv = lognpdf(theta(2),1.6,0.8); %median 5 days, 95% CI [1,24]
    p_alpha = lognpdf(theta(3),0,0.8); %median 1, 95% CI [0.2,5]
    p_beta = lognpdf(theta(4),0.7,0.8); %median 2, 95% CI [0.4,10]
    
    p = p_p_E*p_mu_inv*p_alpha*p_beta;
end