%PLOT_KDE_TO_DATA.m

%% KDE Model Selection and Optimal Bandwidth Estimation

clear; clc; close all;

% --- 0. Load Data Here ---


N=length(Data);

% --- 1. Optimal Bandwidth Estimation (h) ---
% The optimal bandwidth (h) for KDE is chosen by minimizing the AMISE.
% The Normal Reference Rule (a standard, simple method) for h:
% h_optimal = sigma * (4 / (3 * N))^(1/5)
% where sigma is an estimate of the standard deviation.

sigma_data = std(Data);

% Use a robust estimate of sigma, often the minimum of the sample standard 
% deviation and the normalized Interquartile Range (IQR), as recommended by Scott/Silverman.
IQR = iqr(Data);
sigma_robust = min(sigma_data, IQR / 1.349);

% AMISE Optimal Bandwidth formula (Normal Reference Rule)
h_AMISE = sigma_robust * (4 / (3 * N))^(1/5);

fprintf('\nKDE Bandwidth Estimation:\n');
fprintf('  Robust Sigma Estimate: %.4f\n', sigma_robust);
fprintf('  Optimal Bandwidth (h_AMISE): %.4f\n', h_AMISE);

% --- 2. Kernel Density Estimation (KDE) ---
% MATLAB's 'ksdensity' function performs the KDE. 
% We will compare the default bandwidth (often Scott's rule) with AMISE.

% Define the range for plotting
x_plot = linspace(min(Data) - 0.5, max(Data) + 0.5, 200);

% KDE 1: Using the calculated AMISE optimal bandwidth
[f_amise, x_amise] = ksdensity(Data, x_plot, 'Bandwidth', h_AMISE);

% KDE 2: Using the default MATLAB bandwidth (for comparison)
[f_default, x_default] = ksdensity(Data, x_plot); 
h_default = median(f_default) * 10; % Not the exact h, but used for display.

% --- 3. Parametric Fit (Gamma) for Comparison ---
% Reuse the Gamma MLE from the previous module for comparison
mu_data = mean(Data); 
var_data = var(Data);
alpha_hat = mu_data^2 / var_data;
beta_hat = var_data / mu_data;
f_gamma = gampdf(x_plot, alpha_hat, beta_hat);

% --- 4. Visualization ---

figure;
histogram(Data, 'Normalization', 'pdf', 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'k');
hold on;

% Plot 1: True Gamma PDF
f_true = gampdf(x_plot, alpha_true, beta_true);
plot(x_plot, f_true, 'k--', 'LineWidth', 2, 'DisplayName', 'True Gamma PDF');

% Plot 2: Parametric MLE Fit
plot(x_plot, f_gamma, 'r-', 'LineWidth', 2, 'DisplayName', 'Parametric Gamma MLE');

% Plot 3: KDE with AMISE Optimal Bandwidth
plot(x_amise, f_amise, 'b-', 'LineWidth', 2, 'DisplayName', ['KDE (Optimal h=', num2str(h_AMISE, '%.3f'), ')']);

% Plot 4: KDE with Default Bandwidth
plot(x_default, f_default, 'm:', 'LineWidth', 2, 'DisplayName', 'KDE (Default h)');


title('KDE vs. Parametric Fit for Gamma-Simulated Data');
xlabel('Data Value (x)');
ylabel('Density f(x)');
legend('Location', 'best');
grid on;

fprintf('\nVisualization Complete. Observe the bias/variance trade-off:\n');
fprintf('  - Parametric fit (Gamma) has high bias but low variance.\n');
fprintf('  - KDE (Optimal h) minimizes variance while maintaining flexibility.\n');