function LL = maximum_likelihood(beta, DATA, MODEL, Dis)
% maximum_likelihood - Log Likelihood estimation for variogram models in geostatistics.
%
% This function computes the log-likelihood for a given set of spatial data
% and variogram model parameters. It is primarily used for parameter estimation
% in geostatistical modeling, leveraging various covariance models.
%
% INPUTS:
%   beta  - Vector of model parameters [variance, range, nugget]
%           variance: Controls the sill or the total variance captured by the model.
%           range: Influences how quickly the correlation decays with distance.
%           nugget: Represents microscale variation or measurement error.
%   DATA  - Matrix containing spatial coordinates and associated fluctuations or data values
%           [X, Y, Fluctuations], where:
%           X: X-coordinates of data points.
%           Y: Y-coordinates of data points.
%           Fluctuations: Observed fluctuations or data values at the coordinates.
%   MODEL - Integer specifying the type of covariance model to use for variogram:
%           1. Exponential
%           2. Gaussian
%           3. Spherical
%           ... (extendable to other models as needed)
%   Dis   - Distance metric used for calculating pairwise distances:
%           1: Manhattan distance (L1 norm).
%           2: Euclidean distance (L2 norm, default).
%
% OUTPUT:
%   LL    - Computed log-likelihood value, used for model fitting and parameter estimation.
%
% EXAMPLE USAGE:
%   beta = [1, 2, 0.1];  % Example parameters: [variance, range, nugget]
%   DATA = [X, Y, Fluctuations];  % Data matrix with coordinates and values
%   MODEL = 1;  % Using Exponential model
%   Dis = 2;  % Applying Euclidean distance
%   LL = maximum_likelihood(beta, DATA, MODEL, Dis);  % Calculate log-likelihood

% Extract coordinates and fluctuations from input data matrix
X = DATA(:, 1);  % X-coordinates
Y = DATA(:, 2);  % Y-coordinates
F = DATA(:, 3);  % Fluctuations

n = length(X);  % Number of data points

% Compute the covariance matrix using the specified model and parameters
C = covmat(X, Y, MODEL, beta, Dis);

% Calculate the determinant of the covariance matrix, handling potential issues with zero determinant
if det(C) ~= 0
    DC = log(det(C));   % Compute the log of the determinant if it is non-zero
else
    % If the determinant is zero, fall back on using eigenvalues for numerical stability
    AA = log(eig(C));
    DC = sum(AA); % Sum of the log of eigenvalues as an alternative
end

% Compute the log-likelihood using the formula for multivariate normal distributions
LL = 0.5 * (transpose(F) / C * F) + 0.5 * DC;

end

