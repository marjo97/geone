function [EKT, Xk, Yk, SF, PLG] = ordkrig(Xp, Yp, X, Y, Diak, Param, MODEL, maxdis, minneigh, dist_model, N)
% ordkrig - Ordinary Kriging interpolation for 2D spatial data
%
% This function performs Ordinary Kriging on a set of spatial data points
% and estimates values on a defined grid based on a given variogram model.
% It calculates the kriging weights and provides both the estimated values
% and the kriging standard deviations for each grid point.
%
% INPUTS:
%   Xp, Yp    - Grid coordinates for prediction points (reshaped into vectors)
%   X, Y      - Coordinates of the original data points
%   Diak      - Data values corresponding to the coordinates
%   Param     - Parameters for the variogram model
%   MODEL     - Variogram model type:
%               1. Exponential
%               2. Gaussian
%               3. Spherical
%               4. Power
%               5. Spatial covariance (specific to model 46a)
%               6. SSRF covariance model
%   maxdis    - Maximum distance within which to search for neighbors
%   minneigh  - Minimum number of neighbors required to perform kriging
%   dist_model- Distance metric for variogram calculation (unused here)
%   N         - Distance exponent (e.g., 2 for Euclidean distance)
%
% OUTPUTS:
%   EKT       - Matrix of estimated values at the grid points
%   Xk, Yk    - Unique grid coordinates (after reshaping)
%   SF        - Matrix of kriging standard deviations at the grid points
%   PLG       - Matrix indicating the number of neighbors used for each grid point
%
% NOTES:
%   - The function utilizes various variogram models to calculate covariances.
%   - Only grid points with at least the minimum required number of neighbors are processed.
%   - Kriging weights are computed using the covariance matrices.

% Initialize variables
n = length(X);  % Number of original data points
Xk = unique(Xp);  % Unique X-coordinates of grid
Yk = unique(Yp);  % Unique Y-coordinates of grid
nx = length(Xk);  % Number of unique X grid points
ny = length(Yk);  % Number of unique Y grid points
EKT = NaN(ny, nx);  % Initialize matrix for estimated values
SF = NaN(ny, nx);  % Initialize matrix for kriging standard deviations
PLG = zeros(ny, nx);  % Initialize matrix for neighbor counts

% Variogram model functions
gexp = @(expo, r) expo(3) + expo(1) * (1 - exp(-3 * (r / expo(2))));
ggau = @(gau, r) gau(3) + gau(1) * (1 - exp(-3 * (r / gau(2)).^2));
gpow = @(pow, r) pow(3) + pow(1) * r.^pow(2);
gsphe = @(gsph, r) gsph(3) + gsph(1) * (1 - (1.5 * r ./ gsph(2) - 0.5 * (r ./ gsph(2)).^3) .* (r <= gsph(2)));
gssrf = @(ssrf, r) ssrf(4) + ssrf(1)*(1 - exp(-3*(r/ssrf(2)).^ssrf(3)));

% Compute the covariance matrix for the original data points
distan = sqrt(bsxfun(@minus, X, X').^2 + bsxfun(@minus, Y, Y').^2); % Pairwise distance matrix
switch MODEL
    case 1
        CXol = gexp(Param, distan);
    case 2
        CXol = ggau(Param, distan);
    case 3
        CXol = gsphe(Param, distan);
    case 4
        CXol = gpow(Param, distan);
    case 5
        CXol = spatial_cov(Param, distan);
    case 6
        CXol = gssrf(Param, distan);
end

% Iterate over each grid point for kriging estimation
for LOOPi = 1:ny
    for LOOPj = 1:nx
        xshm = Xk(LOOPj);
        yshm = Yk(LOOPi);

        % Compute distances from the current grid point to all original points
        AP = sqrt((xshm - X).^2 + (yshm - Y).^2);

        % Identify neighbors within the specified distance
        neighbors = find(AP < maxdis & AP > 1000 * eps);
        sgt = length(neighbors);
        PLG(LOOPi, LOOPj) = sgt;

        % Proceed with kriging if enough neighbors are found
        if sgt >= minneigh
            Xg = [X(neighbors), Y(neighbors), Diak(neighbors)];
            Cx = CXol(neighbors, neighbors);
            Cx(:, sgt + 1) = 1;
            Cx(sgt + 1, :) = 1;
            Cx(sgt + 1, sgt + 1) = 0;

            Cu = zeros(sgt + 1, 1); % Initialize cross-covariance vector
            distan = sqrt((xshm - X(neighbors)).^2 + (yshm - Y(neighbors)).^2);

            % Calculate cross-covariance based on model
            switch MODEL
                case 1
                    Cu(1:sgt, 1) = gexp(Param, distan);
                case 2
                    Cu(1:sgt, 1) = ggau(Param, distan);
                case 3
                    Cu(1:sgt, 1) = gsphe(Param, distan);
                case 4
                    Cu(1:sgt, 1) = gpow(Param, distan);
                case 5
                    Cu(1:sgt, 1) = spatial_cov(Param, distan);
                case 6
                    Cu(1:sgt, 1) = gssrf(Param, distan);
            end
            Cu(sgt + 1) = 1;

            % Compute kriging weights (Lagrange multipliers included)
            L = Cx \ Cu;

            % Estimate the value at the grid point using kriging weights
            EKT(LOOPi, LOOPj) = Xg(:, 3)' * L(1:sgt);

            % Calculate the kriging standard deviation
            SF(LOOPi, LOOPj) = sqrt(L(1:sgt)' * Cu(1:sgt) + L(sgt + 1));
        end
    end
end
end

