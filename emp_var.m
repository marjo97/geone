function [Variogram, distance, Pairsperclass] = emp_var(Lags, f, tole, X2, Y2, Diak, flag, MDIS, dist_model, N)
% emp_var - Computes the empirical semivariogram in 2D.
%
% This function calculates the empirical semivariogram for a set of spatial data points in 2D.
% It computes pairwise distances, groups them into lag classes, and calculates semivariance for each class.
% Only lag classes with a sufficient number of pairs (at least 30) are considered for reliable semivariance estimates.
%
% INPUTS:
%   Lags       - Number of lag classes to divide the distances into.
%   f          - Direction in radians (0 to pi). Defines the preferred direction for the variogram.
%   tole       - Tolerance for direction in radians (e.g., pi for isotropic analysis).
%   X2, Y2     - Vectors of X and Y coordinates for the data points.
%   Diak       - Vector of data values corresponding to the coordinates (e.g., fluctuations or observations).
%   flag       - Boolean flag indicating whether to plot the semivariogram (1 = yes, 0 = no).
%   MDIS       - Maximum variogram cutoff distance relative to the maximum pair distance (e.g., 0.5 sets it to half the max distance). Default is 0.425.
%   dist_model - Distance metric for pdist2 function (not used directly here, but kept for compatibility).
%   N          - Exponent of distance, specified as a positive scalar (e.g., 2 for Euclidean distance).
%
% OUTPUTS:
%   Variogram     - Vector of computed semivariance values for each lag class.
%   distance      - Vector of the center distances of each lag class.
%   Pairsperclass - Vector of the number of point pairs in each lag class.
%
% NOTES:
%   - The function is optimized for distances in the range (0, pi/2) due to the directional tolerance setting.
%   - It may not perform well for very small distances due to numerical precision limitations.

% Set default maximum variogram cutoff distance if not provided
if nargin < 8
    MDIS = 0.425;
end

% Ensure MDIS is within valid range
if MDIS > 1
    MDIS = 1;
end

% Convert direction and tolerance from radians to degrees for visualization purposes
f2 = f * 360 / (2 * pi);
f2 = f * 360 / (2 * pi);
f3 = (f - pi/2) * 360 / (2 * pi);
tole2 = tole * 360 / (2 * pi);
Dir1 = num2str(f2);
Dir2 = num2str(f3);
Direct = strcat(' Direction:', Dir1);
Direct2 = strcat(' Direction:', Dir2);
Tole1 = num2str(tole2);
Toleran = strcat(' Tolerance:', Tole1);
Titl = strcat(Direct, ',', Toleran);
Tit2 = strcat(Direct2, ',', Toleran);

% Start timing the function execution to assess performance
t0 = cputime;

% Initialize variables
n = length(Diak);
X = [zeros(length(Diak), 1), X2, Y2, Diak];
[n, ~] = size(X);
nx = X(:, 2);
ny = X(:, 3);
value = X(:, 4);

% Calculate angle tolerance range for directional variogram computation
minf = f - tole;
maxf = f + tole;

% Calculate pairwise distances manually using Minkowski distance
AP = zeros(n, n);
for i = 1:n
    for j = 1:n
        AP(i, j) = sum(abs([X2(i), Y2(i)] - [X2(j), Y2(j)]).^N)^(1/N);
    end
end

% Setup for lag distance calculation
rcl = zeros((Lags + 1), 1);
MAP = MDIS * max(max(AP));
lagd = MAP / Lags / 2;  % Radius of each class
rcl(1) = 0;

% Compute lag class centers
for i = 2:Lags+1
    rcl(i) = rcl(i-1) + 2 * lagd;
end

% Initialize variables for semivariogram calculation
class = zeros(n, n);
t = zeros(n, n);
prwtovar = zeros(2 * (Lags + 1), 1);
pointsinclass = zeros(Lags + 1, 1);

% Calculate semivariogram values
for i = 1:n
    for j = 1:n
        class(i, j) = round(AP(i, j) / (2 * lagd)) + 1;
        if class(i, j) < Lags + 2
            dx = nx(j) - nx(i);
            dy = ny(j) - ny(i) + 0.000001;  % Add small number to avoid division by zero
            t(i, j) = acot(dx / dy + 0.000001);
            if (minf < t(i, j) && t(i, j) < maxf) || (minf - pi < t(i, j) && t(i, j) < maxf - pi)
                if AP(i, j) > eps
                    pointsinclass(class(i, j)) = pointsinclass(class(i, j)) + 1;
                    prwtovar(class(i, j)) = prwtovar(class(i, j)) + (abs(value(i) - value(j)))^0.5;
                end
            end
        end
    end
end

% Variogram calculation using a robust estimator
var = zeros(Lags + 1, 1);
secvar = zeros(Lags + 1, 1);
for i = 1:(Lags + 1)
    secvar(i) = (prwtovar(i) / pointsinclass(i))^4;
    var(i) = secvar(i) / (0.494 / pointsinclass(i) + 0.457) / 2;
    if pointsinclass(i) < 30
        var(i) = NaN;
    end
end

% Calculate lag distances for each class center
distance = zeros(Lags + 1, 1);
for i = 2:Lags + 1
    distance(i) = (i - 1) * 2 * lagd;
end

Variogram = var;
Variogram(1) = NaN;
Pairsperclass = pointsinclass;

% Plot the variogram if the flag is set
if nargin > 6 && flag == 1
    figure;
    axes1 = axes('Parent', gcf, 'fontsize', 16);
    hold(axes1, 'all');
    plot(distance, Variogram, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
    hold on;
    plot(0, 0, 'k.');
    xlabel('distance');
    ylabel('\gamma(r)');
    hold off;
end

% Remove NaN values from output variables for clean results
G = isnan(Variogram);
Variogram(G == 1) = [];
distance(G == 1) = [];
Pairsperclass(G == 1) = [];

% Display computation time
Tol = cputime - t0;
fprintf('Computation time: %.2f seconds\n', Tol);

end

