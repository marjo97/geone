%{
  Script for Geostatistical Analysis and Kriging in Octave
  --------------------------------------------------------
  Author: [Koltsidopoulou Maria-Despoina]
  Version: 8.3.0
  Date: [June 2024]
  Description:
  This script performs various geostatistical analyses, including data loading,
  visualization, parameter estimation using maximum likelihood, variogram modeling,
  cross-validation, error analysis, and kriging interpolation. The script is compatible
  with Octave and includes functions adapted for Octave's syntax and capabilities.

  Requirements:
  - Octave version >= 4.0
  - Necessary toolboxes for statistics and optimization

  Usage:
  This script is intended for educational and research purposes to demonstrate
  geostatistical methods in Octave.

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear workspace, close all figures, and clear command window
clear all;
close all;
clc;

% 1: Exponential
% 5: HCE model
% 6: SSRF model
MODEL = 1; % Default model setting

% Load data from .mat files
base_path = 'NormalizedData2';
data = load(fullfile(base_path, 'DED_GAL.mat'));
LUT = load(fullfile(base_path, 'LUTable.mat')).LUTAL;

% Access specific fields from the loaded data
DED_data = data.DEDAL;
GA_value = data.GA;

% Define distance model and exponent
dist_model = 'minkowski';

% N = 1: manhattan distance
% N = 2: euclidean distance
% N = 15: chebyshev distance
N = 1; % Default exponent setting

% Extract coordinates and values from DED_data
X = DED_data(:, 1); % X-coordinates
Y = DED_data(:, 2); % Y-coordinates
VAL = DED_data(:, 3); % Value column
OVAL = DED_data(:, 4); % original values

% Create a meshgrid for interpolation based on the given step size
step_size = 0.1; % km
[XI, YI] = meshgrid(min(X):step_size:max(X), min(Y):step_size:max(Y));
[gr, st] = size(XI);

% Reshape the meshgrids for easier manipulation
Xp = reshape(XI, gr * st, 1); % Reshaped X coordinates
Yp = reshape(YI, gr * st, 1); % Reshaped Y coordinates

% Get unique coordinates for the grid
Xk = unique(Xp);
Yk = unique(Yp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Visualization and Statistical Analysis
% -------------------------------------------
% This section visualizes the spatial distribution of data points and explores
% the statistical properties of the dataset. It includes:
% - A scatter plot of the data point locations.
% - A histogram and normal probability plot to understand data distribution.
% - Calculation of basic statistical metrics such as min, max, mean, median,
%   standard deviation, skewness, and kurtosis.
% - Empirical cumulative distribution function (CDF) plot compared to a Gaussian CDF.

% Plotting a scatter plot of locations
figure('Name', 'Scatter Plot of Locations', 'Position', [70, 70, 700, 700], 'Color', 'w');
scatter(X, Y, 'x', 'MarkerEdgeColor', 'b', 'LineWidth', 1.5);
axis equal;
grid on;
box on;
xlabel('X (km)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Y (km)', 'FontSize', 12, 'FontWeight', 'bold');
title('Location of Data Points', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'LineWidth', 1.5);
legend('Data Points', 'Location', 'northeast');

% Histogram and Normal Probability Plot
figure('Name', 'Data Distribution', 'Position', [800, 70, 700, 700]);
subplot(2,1,1);
hist(VAL, 20); % Again, using 'hist'
title('Histogram');
xlabel('Values');
ylabel('Frequency');
subplot(2,1,2);
% Manual implementation of a normal probability plot
sorted_VAL = sort(VAL);
p = (1:length(sorted_VAL)) / (length(sorted_VAL) + 1);
theoretical_quantiles = norminv(p, mean(VAL), std(VAL));
plot(theoretical_quantiles, sorted_VAL, 'o');
hold on;
plot([min(theoretical_quantiles), max(theoretical_quantiles)], ...
     [min(sorted_VAL), max(sorted_VAL)], 'r-', 'LineWidth', 2);
hold off;
xlabel('Theoretical Quantiles');
ylabel('Sample Quantiles');
title('Normal Probability Plot');
grid on;

% Data Statistics
% Calculate basic statistics for the dataset
stats = struct('min', min(VAL), 'max', max(VAL), 'mean', mean(VAL), ...
               'median', median(VAL), 'std', std(VAL), 'skewness', skewness(VAL), ...
               'kurtosis', kurtosis(VAL));
disp(stats);

% Empirical CDF and comparison with Gaussian CDF
% Plotting empirical CDF and comparing with Gaussian CDF
figure('Name', 'Empirical and Gaussian CDF', 'Position', [70, 70, 700, 700]);
[F, z] = ecdf(VAL);  % Empirical CDF
stairs(z, F, 'k', 'LineWidth', 2);
hold on;
y = normcdf(z, stats.mean, stats.std);  % Gaussian CDF
plot(z, y, 'r', 'LineWidth', 1);
hold off;
xlabel('Value');
ylabel('Cumulative Probability');
legend('Empirical', 'Gaussian', 'Location', 'southeast');
title('CDF Comparison');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Maximum Likelihood Estimation Setup
% -----------------------------------
% This section prepares for parameter estimation using the Maximum Likelihood
% method. The initial parameter estimates, lower and upper bounds for the
% parameters, and optimization options are defined. Depending on the chosen
% model, specific initial estimates and bounds are set. The objective function
% is wrapped to include additional data needed for the likelihood calculation.
% Finally, constrained optimization is performed using Sequential Quadratic
% Programming (SQP) to estimate the parameters that maximize the likelihood.

% Initial parameter estimates and bounds
beta0 = [0.9, 0.1, eps];
LB = [0.5, 0.05, eps];  % Lower bounds
UB = [2, 0.5, 0.5];  % Upper bounds

% Adjust parameters if a different model is selected
if MODEL == 6
    beta0 = [14.5, 2, 1.195, eps];
    LB = [0.5, 1.8, 0.05/3, eps];
    UB = [20, 2.15, 5/3, 0.5];
end

% Optimization options
options = optimset('display','final','TolX',1e-4,'TolFun',1e-3,...
   'MaxFunEvals',1e5,'MaxIter',1e5);

% Define a wrapper function to pass additional parameters to LLik3
objective_function = @(beta) maximum_likelihood(beta, [X Y VAL], MODEL, N);

% Perform parameter estimation using constrained optimization with sqp
[Param, fval, info, output] = sqp(beta0, objective_function, [], [], LB, UB);
disp('Optimized Parameters:');
disp(Param);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model fitting based on variogram
% --------------------------------
% This section fits different variogram models to the empirical semivariogram
% data computed from the input spatial data. Depending on the selected model,
% the script calculates the theoretical variogram, compares it with the empirical
% variogram, and plots the results for visual analysis.
% The error (ERR) between the modeled and empirical variogram is also computed
% to assess the goodness of fit for each model.
Nc = 15;
Nugesti = 1;  Indicator for including a nugget effect
Param2 = Param; Param2(3) = Param(3); % Ensure consistency in parameters

% Compute the empirical semivariogram
[varioz, lags, pairs] = emp_var(15, pi/2, pi/2, X, Y, VAL, 0, 0.5, dist_model, N);

% Switch-case to handle different variogram models
switch MODEL
    case 1
        % DExponential model
        modexpon = @(betaexp1, x) betaexp1(3) + betaexp1(1) * (1 - exp(-x / betaexp1(2)));
        lagsn = linspace(0, max(lags), 200); % Generate lag values for plotting
        % Calculate the theoretical variogram using the exponential model
        variom = modexpon(Param, lagsn);
        variom_lc = modexpon(Param, lags);
        % Compute error between empirical and theoretical variogram
        ERR = sum((variom_lc - varioz).^2);
        % Plot the empirical and modeled variogram
        fig = figure;
        axes1 = axes('Parent', fig, 'FontSize', 16, 'FontWeight', 'bold', 'Box', 'on');
        hold(axes1, 'all');
        grid(axes1, 'on');
        plot(lags, varioz, 'o', 'MarkerSize', 9, 'MarkerEdgeColor', 'b', 'LineWidth', 1.5);
        plot(lagsn, variom, 'r', 'LineWidth', 2);
        plot(eps, eps, '.', 'MarkerSize', 12, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
        xlabel('Distance r', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
        ylabel('\gamma(r)', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
        legend({'Empirical', 'Exponential Model'}, 'Location', 'northeast', 'FontSize', 14, 'TextColor', 'k', 'Box', 'off');
        set(axes1, 'FontSize', 16, 'FontWeight', 'bold', 'LineWidth', 1.5, 'GridColor', [0.8, 0.8, 0.8], 'GridAlpha', 0.7);
        axis tight;
        hold off;
    case 5
        % HCE covariance model
        modhce = @(beta, x) spatial_cov(beta, eps) - spatial_cov(beta, x);
        lagsn = linspace(eps*1000, max(lags), 200); % Generate lag values for plotting
        % Calculate the theoretical variogram using the HCE model
        [variom] = modhce(Param, lagsn);
        [variom_lc] = modhce(Param, lags);
        % Compute error between empirical and theoretical variogram
        ERR = sum((variom_lc - varioz).^2);
        % Plot the empirical and modeled variogram
        fig = figure;
        axes1 = axes('Parent', fig, 'FontSize', 16, 'FontWeight', 'bold', 'Box', 'on');
        hold(axes1, 'all');
        grid(axes1, 'on');
        plot(lags, varioz, 'o', 'MarkerSize', 9, 'MarkerEdgeColor', 'b', 'LineWidth', 1.5);
        plot(lagsn, variom, 'r', 'LineWidth', 2);
        plot(eps, eps, '.', 'MarkerSize', 12, 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g');
        xlabel('Distance r', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
        ylabel('\gamma(r)', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
        legend({'Empirical', 'HCE Covariance Model'}, 'Location', 'northeast', 'FontSize', 14, 'TextColor', 'k', 'Box', 'off');
        set(axes1, 'FontSize', 16, 'FontWeight', 'bold', 'LineWidth', 1.5, 'GridColor', [0.8, 0.8, 0.8], 'GridAlpha', 0.7);
        axis tight;
        hold off;
    case 6
        % SSRF covariance model
        lagsn = linspace(eps*1010, max(lags), 200); % Generate lag values for the variogram model plot
        % Calculate SSRF covariance values
        [~, varsrf] = SSRF_cov(3, eps*1000, Param, 0);
        [~, variosrf] = SSRF_cov(3, lagsn, Param, 0);
        % Compute the model variogram from SSRF covariances
        variom = varsrf - variosrf;
        if Nugesti == 1
            variom = variom + Param2(4);
        end
        % Plot the empirical and modeled variogram
        figure;
        hold on;
        plot(lags, varioz, 'o', 'MarkerSize', 9, 'MarkerEdgeColor', 'b', 'LineWidth', 1.5);
        plot(lagsn, variom, 'r', 'LineWidth', 2);
        xlabel('Distance r', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
        ylabel('\gamma(r)', 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
        legend({'Empirical', 'SSRF Covariance Model'}, 'Location', 'northeast', 'FontSize', 14, 'TextColor', 'k', 'Box', 'off');
        grid on;
        set(gca, 'GridColor', [0.8, 0.8, 0.8], 'GridAlpha', 0.7);
        set(gca, 'FontSize', 16, 'FontWeight', 'bold', 'LineWidth', 1.5, 'Box', 'on');
        axis tight;
        hold off;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cross-validation and error analysis
% -----------------------------------
% Perform cross-validation to evaluate the model's predictive accuracy.
% Identify and remove outliers based on cross-validation errors for
% improved analysis. Compute error metrics after back-transforming
% cross-validated estimates to the original data scale, including Mean
% Absolute Error (MAE), Root Mean Squared Error (RMSE), correlation
% coefficient, and Mean Error (ME). Display these metrics for assessing
% the model's performance and reliability.
NERAD = 70 / 100; % 70% of range
[ECV, ~, PLG, Err, LTot] = crossval(X, Y, VAL, Param, MODEL, NERAD, [], dist_model, N);


% Detecting outliers based on cross-validation errors
KATS = abs(ECV - VAL) > 5;
if sum(KATS) >= 5
    disp('*****************')
    disp('Many OUTLIERS! ! ! !')
    disp('*****************')
    OUTLIERS = sum(KATS);
end

% Remove outliers for further analysis
ECV2 = ECV;
VAL2 = VAL;
OVAL2 = OVAL;
ECV2(KATS) = NaN;
G = isnan(ECV2);
ECV2(G) = [];
VAL2(G) = [];
OVAL2(G) = [];

% Back-transformed error metrics calculation and display
Ka = LUT(:,2);  % Knnsearch output array
Or = LUT(:,1);  % Original values array from LUT
D = knnsearch(Ka, ECV2);  % Find nearest neighbors
BTECV = Or(D);  % Backtransformed Estimated Cumulative Volume
ErrBackTransformed = BTECV - OVAL2;  % Errors after back-transformation
MAE_BackTransformed = mean(abs(ErrBackTransformed));  % Mean Absolute Error for back-transformed data
RMSE_BackTransformed = sqrt(mean(ErrBackTransformed.^2));  % Root Mean Squared Error for back-transformed data
rho_BackTransformed = corr(BTECV, OVAL2);  % Correlation coefficient for back-transformed data
ME_BackTransformed = mean(ErrBackTransformed);  % Mean Error for back-transformed data

% Display back-transformed results
disp('Back Transformed Results:');
fprintf('Back Transformed MAE: %f\n', MAE_BackTransformed);
fprintf('Back Transformed RMSE: %f\n', RMSE_BackTransformed);
fprintf('Back Transformed Correlation Coefficient: %f\n', rho_BackTransformed);
fprintf('Back Transformed Mean Error: %f\n', ME_BackTransformed);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform Kriging to Interpolate Spatial Data
% -------------------------------------------------
% This section performs ordinary kriging interpolation on the data
% points to predict values at specified grid locations.
NERAD = 100 / 100;
[EKT, Xk, Yk, ~, PLG] = ordkrig(Xp, Yp, X, Y, VAL, Param, MODEL, NERAD, 3, dist_model, N);

[m, n] = size(EKT);
GG = isnan(EKT);
EKT = reshape(EKT, m * n, 1);

D2 = knnsearch(Ka, EKT);
K = Or(D2);
K(K < 1) = 1;
B_EKT = log(K);
B_EKT = reshape(B_EKT, m, n);
B_EKT(GG) = nan;

% Plot the Kriging results using the modified function
kriging_after_bt(Xk, Yk, B_EKT);
