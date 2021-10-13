function p = prior_fun_indep(theta)
    
    % Calculate the prior density for the independent transmission and
    % symptoms model with estimated parameters theta
    
    p_mean = lognpdf(theta(1),1.6,0.35); %median 5, 95% CI [2.5,10]
    p_sd = lognpdf(theta(2),0.7,0.65); %median 2, 95% CI [0.6,7]
    p_beta = lognpdf(theta(3),0.7,0.8); %median 2, 95% CI [0.4,10]
    
    p = p_mean*p_sd*p_beta;
end