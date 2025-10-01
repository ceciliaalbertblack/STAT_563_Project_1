function [AIC, SBC, ICOMP_F, theta_hat, LogL_max] = fit_and_select_pdf(Data, logL_func, theta_0, dist_name)
% FIT_AND_SELECT_PDF
%   Performs MLE, calculates Information Criteria (AIC, SBC, ICOMP_F), 
%   and finds the estimated asymptotic covariance matrix (Sigma_hat)
%   using the Hessian from FMINUNC.
%
% Inputs:
%   Data      : The observed data vector.
%   logL_func : An anonymous function handle for the NEGATIVE Log-Likelihood.
%               Syntax: @(theta) -sum(log(pdf_func(Data, theta(1), ...)))
%   theta_0   : Initial guess vector for the parameters.
%   dist_name : String name of the distribution (for display).
%
% Outputs:
%   AIC, SBC, ICOMP_F, theta_hat, LogL_max

    N = length(Data);
    m = length(theta_0);

    % --- 1. Maximum Likelihood Estimation (MLE) using FMINUNC ---
    % FMINUNC is used because it can return the Hessian matrix H, needed for ICOMP_F.
    % H at the MLE is approximately -Fisher Information Matrix (I).
    try
        options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'off');
        [theta_hat, negLogL_min, ~, ~, ~, H_hat] = fminunc(logL_func, theta_0, options);
        
        LogL_max = -negLogL_min;
        
        fprintf('\n--- Results for %s ---\n', dist_name);
        fprintf('  Parameters (m): %d\n', m);
        fprintf('  Theta_hat: [%s]\n', num2str(theta_hat, '%.4f '));
        fprintf('  Log-Likelihood (L_max): %.4f\n', LogL_max);
        
    catch ME
        fprintf('\n--- ERROR: %s ---\n', dist_name);
        fprintf('  Optimization failed (requires Optimization Toolbox for fminunc).\n');
        fprintf('  Error: %s\n', ME.message);
        AIC = NaN; SBC = NaN; ICOMP_F = NaN; theta_hat = NaN; LogL_max = NaN;
        return;
    end

    % --- 2. Calculate AIC and SBC ---
    
    % AIC = -2*L_max + 2*m
    AIC = -2 * LogL_max + 2 * m;
    
    % SBC = -2*L_max + m*ln(N)
    SBC = -2 * LogL_max + m * log(N);
    
    % --- 3. Calculate ICOMP_F ---
    
    % Estimated Asymptotic Covariance Matrix: Sigma_hat = I^-1 ~ -H^-1
    % Ensure the Hessian is positive definite (Hessian of NLL should be P.D.)
    try
        Sigma_hat = inv(H_hat); 
    catch
        % If H_hat is singular (often means model is over-specified), use pseudoinverse
        Sigma_hat = pinv(H_hat);
    end
    
    % Ensure Sigma_hat is symmetric and deal with potential numerical issues
    Sigma_hat = (Sigma_hat + Sigma_hat') / 2;
    
    % Calculate eigenvalues and mean eigenvalue
    eigen_Sigma = eig(Sigma_hat);
    
    % Filter out non-real or near-zero eigenvalues if they are present due to numerical issues
    eigen_Sigma = real(eigen_Sigma(abs(eigen_Sigma) > 1e-10)); 
    
    if isempty(eigen_Sigma) || any(eigen_Sigma < 0)
        warning('ICOMP: Negative or missing eigenvalues encountered. ICOMP_F set to NaN.');
        ICOMP_F = NaN;
        
    else
        m_eff = length(eigen_Sigma); % Effective number of parameters for ICOMP
        mean_eigen = mean(eigen_Sigma);
        
        % C1F Penalty Formula:
        C1F_penalty = (1 / (4 * mean_eigen^2)) * sum((eigen_Sigma - mean_eigen).^2);
        
        % ICOMP_F = -2*L_max + C1F_penalty
        ICOMP_F = -2 * LogL_max + C1F_penalty;
    end

    fprintf('  AIC:   %.4f\n', AIC);
    fprintf('  SBC:   %.4f\n', SBC);
    fprintf('  ICOMP_F: %.4f\n', ICOMP_F);
end
