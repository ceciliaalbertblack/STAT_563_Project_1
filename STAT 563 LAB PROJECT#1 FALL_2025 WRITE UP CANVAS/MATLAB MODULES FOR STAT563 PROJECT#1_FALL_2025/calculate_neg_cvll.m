    function nll = calculate_neg_cvll(data, h)
        % Calculates the negative Log-Likelihood using Leave-One-Out Cross-Validation
        N_data = length(data);
        log_likelihood_sum = 0;
        
        % The ksdensity function handles the leave-one-out implicitly when 
        % its 'kernel' option is set to a CV method, but a manual loop is 
        % the most explicit demonstration of the CV principle:
        
        for i = 1:N_data
            % Data set excluding the i-th point
            data_minus_i = data([1:i-1, i+1:end]);
            
            % Estimate density at x_i using the rest of the data
            % KSDENSITY uses a 'normal' kernel by default
            f_hat_at_xi = ksdensity(data_minus_i, data(i), 'Bandwidth', h);
            
            % Avoid log(0)
            if f_hat_at_xi > eps
                log_likelihood_sum = log_likelihood_sum + log(f_hat_at_xi);
            else
                log_likelihood_sum = log_likelihood_sum - 10; % Heavy penalty for near zero
            end
        end
        nll = -log_likelihood_sum; % Return negative CVLL
    end