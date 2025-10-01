%KDE_BANDWIDTH_SELECTION.m

%% MATLAB Script: Hybrid KDE Scoring (CVLL + Normal FIM Penalty)

clear; clc; close all;

% Load or import the data here

% --- 0. Load Data Here ---



N=length(Data);

% --- 1. CVLL Optimization for Optimal h ---

% Define the NEGATIVE CVLL loss function (to be minimized)
negCVLL_loss = @(h) calculate_neg_cvll(Data, h);

% Initial guess for h (AMISE Normal Reference Rule)
std_data = std(Data);
h0 = std_data * (4 / (3 * N))^(1/5); 

% Use fminsearch to find the h that MINIMIZES the negative CVLL
options = optimset('Display', 'off');
[h_opt_cvll, negCVLL_min] = fminsearch(negCVLL_loss, h0, options);

CVLL_max = -negCVLL_min;
m = 2; % Assume 2 parameters: mu and sigma

fprintf('\nKDE Hybrid Scoring Results:\n');
fprintf('  Optimal Bandwidth (h_CVLL): %.4f\n', h_opt_cvll);
fprintf('  CV Log-Likelihood (L_max proxy): %.4f\n', CVLL_max);


% --- 2. Hybrid Information Criteria Calculation ---

% A. AIC (Hybrid)
% AIC = -2*CVLL_max + 2*m
AIC_hybrid = -2 * CVLL_max + 2 * m;

% B. SBC/BIC (Hybrid)
% SBC = -2*CVLL_max + m*ln(N)
SBC_hybrid = -2 * CVLL_max + m * log(N);

% C. ICOMP_F (Hybrid - Using Normal FIM Penalty)
mu_hat = mean(Data);
sigma_hat = std_data;

% FIM-Inverse (Asymptotic Covariance Matrix) for Normal Distribution (N observations)
% Sigma_hat = 1/N * [sigma^2, 0; 0, sigma^2/2]
Sigma_hat = (1/N) * [sigma_hat^2, 0; 0, 2*sigma_hat^4];

% Calculate ICOMP_F Penalty (C1F)
eigen_Sigma = eig(Sigma_hat);
mean_eigen = mean(eigen_Sigma);

% C1F Penalty Formula:
C1F_penalty = (1 / (4 * mean_eigen^2)) * sum((eigen_Sigma - mean_eigen).^2);

% ICOMP_F = -2*CVLL_max + C1F_penalty
ICOMP_F_hybrid = -2 * CVLL_max + C1F_penalty;

fprintf('\nHybrid Information Criteria Scores (Lower is Better):\n');
fprintf('  AIC Hybrid:   %.4f\n', AIC_hybrid);
fprintf('  SBC Hybrid:   %.4f\n', SBC_hybrid);
fprintf('  ICOMP_F Hybrid: %.4f\n', ICOMP_F_hybrid);

% --- Visualization ---
figure;
histogram(Data, 'Normalization', 'pdf', 'FaceColor', [0.7 0.9 0.7]);
hold on;
x_plot = linspace(min(Data), max(Data), 200);

% Plot KDE with CVLL-optimized bandwidth
[f_cvll, ~] = ksdensity(Data, x_plot, 'Bandwidth', h_opt_cvll);
plot(x_plot, f_cvll, 'r-', 'LineWidth', 2, 'DisplayName', ['KDE (Hybrid Optimal h=', num2str(h_opt_cvll, '%.4f'), ')']);

title('KDE Optimal Bandwidth Selection using Hybrid IC Scoring');
xlabel('Data Value');
ylabel('Density');
legend('Location', 'best');
grid on;


%writematrix(Data, 'Raw_Project_Data.xlsx');