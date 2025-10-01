%RUN_PROJECT.m

% --------------------------------------------------------------------------
% --- MAIN SCRIPT TO RUN THE FUNCTION AND COMPARE DISTRIBUTIONS ---
% --------------------------------------------------------------------------
% The code below is a separate script to demonstrate calling the function
% multiple times for different candidate distributions.

clear; clc; close all;

% --- 0. Load Data Here  ---


N=length(Data);

% --- 1. Define Log-Likelihoods and Initial Guesses for ALL Candidates ---

% General Data Statistics for Initial Guesses
mu_data = mean(Data); 
var_data = var(Data);
std_data = std(Data);
logData = log(Data);

% --- Candidate 1: Gamma Distribution (Parameters: [alpha, beta]) ---
gamma_logL = @(theta) -sum(log(gampdf(Data, theta(1), theta(2))));
theta0_gamma = [mu_data^2 / var_data, var_data / mu_data]; 

% --- Candidate 2: Log-Normal Distribution (Parameters: [mu, sigma]) ---
logL_lognormal = @(theta) -sum(log(lognpdf(Data, theta(1), theta(2))));
theta0_lognormal = [mean(logData), std(logData)];

% --- Candidate 3: Weibull Distribution (Parameters: [scale, shape]) ---
logL_weibull = @(theta) -sum(log(wblpdf(Data, theta(1), theta(2))));
theta0_weibull = [mu_data, 1]; % Mean and shape=1 (Exponential)

% --- Candidate 4: Normal (Gaussian) Distribution (Parameters: [mu, sigma]) ---
logL_normal = @(theta) -sum(log(normpdf(Data, theta(1), theta(2))));
theta0_normal = [mu_data, std_data];

% --- Candidate 5: Extreme Value (Gumbel) Distribution (Parameters: [mu, sigma]) ---
logL_ev = @(theta) -sum(log(evpdf(Data, theta(1), theta(2))));
theta0_ev = [mu_data, std_data]; % Simple moment guess

% --- Candidate 6: Generalized Extreme Value (GEV) (Parameters: [mu, sigma, xi]) ---
logL_gev = @(theta) -sum(log(gevpdf(Data, theta(1), theta(2), theta(3))));
theta0_gev = [mu_data, std_data, 0.1]; % Add small shape parameter guess

% --- Candidate 7: Generalized Pareto (GP) Distribution (Parameters: [k, sigma, theta]) ---
% Data must be greater than the threshold (theta).
logL_gppdf = @(theta) -sum(log(gppdf(Data, theta(1), theta(2), theta(3))));
theta0_gppdf = [0.1, std_data, min(Data)]; % k_0=0.1, sigma_0=std, theta_0=min(Data)


% --- 2. Fit Models and Collect Results ---
model_names = {'Gamma', 'Log-Normal', 'Weibull', 'Normal', 'Extreme Value', 'GEV', 'GP'};
logL_funcs = {gamma_logL, logL_lognormal, logL_weibull, logL_normal, logL_ev, logL_gev, logL_gppdf};
theta_initials = {theta0_gamma, theta0_lognormal, theta0_weibull, theta0_normal, theta0_ev, theta0_gev, theta0_gppdf};

num_models = length(model_names);
Results = table('Size', [num_models, 5], ...
    'VariableTypes', {'string', 'double', 'double', 'double', 'double'}, ...
    'VariableNames', {'Model', 'AIC', 'SBC', 'ICOMP_F', 'L_max'});

for i = 1:num_models
    [AIC, SBC, ICOMP_F, ~, L_max] = fit_and_select_pdf(Data, logL_funcs{i}, theta_initials{i}, model_names{i});
    
    Results.Model(i) = model_names{i};
    Results.AIC(i) = AIC;
    Results.SBC(i) = SBC;
    Results.ICOMP_F(i) = ICOMP_F;
    Results.L_max(i) = L_max;
end

%% 3. Final Comparison and Visualization

fprintf('\n\n======================================================\n');
fprintf('FINAL MODEL SELECTION SUMMARY (Minimum Value is Best)\n');
fprintf('======================================================\n');
disp(Results);

% Find best model for each criterion
[~, idx_aic] = min(Results.AIC);
[~, idx_sbc] = min(Results.SBC);
[~, idx_icomp] = min(Results.ICOMP_F);

fprintf('\nOptimal Model by AIC:   %s (Value: %.4f)\n', Results.Model{idx_aic}, Results.AIC(idx_aic));
fprintf('Optimal Model by SBC:   %s (Value: %.4f)\n', Results.Model{idx_sbc}, Results.SBC(idx_sbc));
fprintf('Optimal Model by ICOMP: %s (Value: %.4f)\n', Results.Model{idx_icomp}, Results.ICOMP_F(idx_icomp));