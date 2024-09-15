function [Cx] = covmat(X, Y, MODEL, Param, DIS)
% covmat - Computes the covariance matrix for a 2D field using a specified covariance model.
%
% This function calculates the covariance matrix for a set of spatial data points
% based on their coordinates (X, Y) and a chosen covariance model. The function supports
% multiple models commonly used in geostatistics, allowing for flexible modeling of spatial
% relationships. The covariance matrix is crucial for various geostatistical analyses,
% including Kriging and spatial prediction.
%
% INPUTS:
%   X      - Vector of X-coordinates of the data points.
%   Y      - Vector of Y-coordinates of the data points.
%   MODEL  - Integer specifying the covariance model to use:
%            1. Exponential model (default)
%            2. Gaussian model
%            3. Spherical model
%            4. NewMODEL 44a (Parameters are Sill and range, not s^2 and b) [Not Implemented]
%            5. NewMODEL 46a (Custom model based on SpatialCov_b function)
%            6. SSRF model [h0, h1, correlation length, nugget]
%   Param  - Vector of parameters for the selected covariance model:
%            [Stdev, range, Nugget]. Specific to each model.
%   DIS    - Minkowski distance exponent:
%            1 for Manhattan distance
%            2 for Euclidean distance
%
% OUTPUT:
%   Cx - The computed covariance matrix (n x n) for n data points.
%
% EXAMPLE USAGE:
%   X = rand(30, 1) * 3;  % Generate random X-coordinates
%   Y = rand(30, 1) * 3;  % Generate random Y-coordinates
%   Cx = covmat(X, Y, 1, [1, 2.5, 0.2], 2);  % Compute covariance using the exponential model with Euclidean distance
%
% NOTE:
%   - The function requires all inputs to be properly defined and consistent with the model specifications.
%   - Currently, MODEL 4 is not implemented and will result in an error if selected.

% Combine X and Y coordinates into a single matrix for distance calculation
LOC = [X Y];

% Calculate pairwise distances using the Minkowski distance metric
n = size(LOC, 1);
DMAT = zeros(n, n);

% Compute the distance matrix using the specified Minkowski exponent (DIS)
for i = 1:n
    for j = 1:n
        DMAT(i, j) = sum(abs(LOC(i, :) - LOC(j, :)).^DIS)^(1/DIS);
    end
end

% Define inline functions for the covariance models
gexp = @(expo, r) expo(1) * exp(-3 * (r / expo(2))); % Exponential model
ggau = @(gau, r) gau(1) * exp(-3 * (r / gau(2)).^2); % Gaussian model
gsph = @(sph, r) sph(1) * (1 - 1.5 * r / sph(2) + 0.5 * (r / sph(2)).^3) .* (r < sph(2)); % Spherical model

% Initialize the covariance matrix based on the selected model
switch MODEL
    case 1  % Exponential model
        Cx = gexp(Param, DMAT);
    case 2  % Gaussian model
        Cx = ggau(Param, DMAT);
    case 3  % Spherical model
        Cx = gsph(Param, DMAT);
    case 5  % HCE model
        b = Param(2) / 2.506;
        V = 6.281 * b^2 * Param(1);
        Cx = spatial_cov([V, b, 0], DMAT);  % Using nugget = 0
    case 6  % SSRF model
        [~, Cx] = SSRF_cov(3, DMAT, Param, 0);
    otherwise
        error('Model not supported or not implemented.'); % Handle unsupported models
end

% Add the nugget effect to the diagonal of the covariance matrix
Cx = Cx + eye(n) * Param(3);

end

