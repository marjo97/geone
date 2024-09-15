function [EKT, SF, PLG, Err, LTot] = crossval(X, Y, Diak, Param, MODEL, maxdis, ArWei, dist_model, N)
% crossval - Performs cross-validation for kriging models.
%
% This function conducts a leave-one-out cross-validation for kriging
% to evaluate the performance of various variogram models. It calculates
% prediction errors, kriging standard deviations, and tracks the number
% of neighbors used for each prediction.
%
% INPUTS:
%   X, Y      - Vectors of coordinates for the data points.
%   Diak      - Data values corresponding to the coordinates.
%   Param     - Parameters for the variogram model.
%   MODEL     - Variogram model type:
%               1. Exponential
%               2. Gaussian
%               3. Spherical
%               4. Power
%   maxdis    - Maximum distance within which to search for neighbors.
%   ArWei     - Optional weight array for weighted error metrics.
%   dist_model- Distance metric (unused in this function, included for compatibility).
%   N         - Distance exponent (unused in this function, included for compatibility).
%
% OUTPUTS:
%   EKT       - Vector of estimated values at the data points after leave-one-out.
%   SF        - Vector of kriging standard deviations for each point.
%   PLG       - Vector indicating the number of neighbors used for each prediction.
%   Err       - Array of error metrics [ME, MAE, RMSE, MaxAE, Rho].
%   LTot      - Cell array of kriging weights for each data point.
%
% NOTES:
%   - The function performs kriging without a specific distance metric or exponent.
%   - It supports weighted and unweighted error calculations.
%   - The distance calculations are done manually to avoid self-interaction.

% Set default variogram model if not specified
if nargin < 5
    MODEL = 1;
end

LTot = cell(length(X), 1); % Cell array to store kriging weights
mindis = 1000 * eps; % Minimum distance to avoid self-interaction
XT = [X Y Diak];
n = length(X); % Number of data points
EKT = zeros(n, 1);

% Define variogram model functions
gexp = @(expo, r) expo(3) + expo(1) * (1 - exp(-3 * (r / expo(2))));
ggau = @(gau, r) gau(3) + gau(1) * (1 - exp(-3 * (r / gau(2)).^2));
gpow = @(pow, r) pow(3) + pow(1) * r.^pow(2);
gsphe = @(gsph, r) gsph(3) + gsph(1) * (1 - (1.5 * r ./ gsph(2) - 0.5 * (r ./ gsph(2)).^3) .* (r <= gsph(2)));

% Precompute pairwise distances
D = sqrt(bsxfun(@minus, X(:), X(:)').^2 + bsxfun(@minus, Y(:), Y(:)').^2);

% Determine if weights are provided for error metrics
if nargin > 6 && length(ArWei) == length(X)
    Weig = 1;
else
    Weig = 0;
end

% Set default maximum distance if not specified
if nargin < 6 && MODEL ~= 4
    maxdis = 0.5;
elseif nargin < 6 && MODEL == 4
    maxdis = 0.1 * max(max(D));
end

maxdis = 0.3;

% Ensure max distance does not exceed one-third of maximum pairwise distance
if maxdis > 0.33 * max(max(D))
    maxdis = 0.33 * max(max(D));
end

clear D % Clear distance matrix to save memory

% Loop through each data point for leave-one-out cross-validation
for LOOP = 1:n
  % Display progress for every 1000 points
    if mod(LOOP, 1000) == 0
        disp(['Points investigated so far: ', num2str(LOOP)]);
    end

    % Extract coordinates of the current data point
    xshm = XT(LOOP, 1);
    yshm = XT(LOOP, 2);

    % Calculate distances from the current point to all other points
    AP = sqrt((xshm - XT(:, 1)).^2 + (yshm - XT(:, 2)).^2); % Distance between the validated point and all others

    % Initialize variables for neighbors
    sgt = 0;
    x = NaN;
    y = NaN;
    v = NaN;

    % Identify valid neighbors within the specified distance
    for i = 1:n
        if AP(i) > mindis && AP(i) < maxdis
            sgt = sgt + 1;
            x(sgt, 1) = XT(i, 1);
            y(sgt, 1) = XT(i, 2);
            v(sgt, 1) = XT(i, 3);
        end
    end

    PLG(LOOP, 1) = sgt; % Store number of neighbors for current point

    % Only proceed if there are neighbors available
    if sgt > 0
        Cx = zeros(sgt + 1, sgt + 1);
        distan = sqrt(bsxfun(@minus, x, x').^2 + bsxfun(@minus, y, y').^2);
        switch MODEL
            case 1
                Cx(1:sgt, 1:sgt) = gexp(Param, distan);
            case 2
                Cx(1:sgt, 1:sgt) = ggau(Param, distan);
            case 3
                Cx(1:sgt, 1:sgt) = gsphe(Param, distan);
            case 4
                Cx(1:sgt, 1:sgt) = gpow(Param, distan);
        end

        % Add Lagrange multipliers to the covariance matrix
        Cx(:, sgt + 1) = 1;
        Cx(sgt + 1, :) = 1;
        Cx(sgt + 1, sgt + 1) = 0;
    end

    if sgt > 0
        Cu = zeros(sgt + 1, 1);
        % Calculate cross-covariance vector
        distan = sqrt((xshm - x).^2 + (yshm - y).^2);
        switch MODEL
            case 1
                Cu(1:sgt, 1) = gexp(Param, distan);
            case 2
                Cu(1:sgt, 1) = ggau(Param, distan);
            case 3
                Cu(1:sgt, 1) = gsphe(Param, distan);
            case 4
                Cu(1:sgt, 1) = gpow(Param, distan);
        end
        Cu(sgt + 1) = 1;
    end

    % Compute kriging weights
    if sgt > 0
        L = Cx \ Cu;
        LTot{LOOP} = L;
    end

    % Calculate kriging estimate and standard deviation
    EKT(LOOP, 1) = NaN;
    SF(LOOP, 1) = NaN;
    if sgt > 0
        EKT(LOOP, 1) = v' * L(1:sgt);
        SF(LOOP, 1) = sqrt(L(1:sgt)' * Cu(1:sgt) + L(sgt + 1));
    end
end

% Remove NaN values for error calculation
G = isnan(EKT);
Diakk = Diak(G == 0);
EKTk = EKT(G == 0);

% Compute errors based on weighted or unweighted metrics
if Weig == 1
    % Weighted error metrics
    ArWei = ArWei(G == 0);

    Er = (Diakk - EKTk);
    SqEr = (Diakk - EKTk).^2;

    disp('Weighted Mean Error');
    ME = Er' * ArWei;

    disp('Weighted Mean Absolute Error');
    MAE = (abs(Er)') * ArWei;

    disp('Weighted Root Mean Square Error');
    RMSE = sqrt(SqEr' * ArWei);

    disp('Weighted correlation coefficient');
    A = sum(Diakk .* EKTk .* ArWei) - sum(EKTk .* ArWei) * sum(Diakk .* ArWei);
    B = sqrt(sum(Diakk.^2 .* ArWei) - sum(Diakk .* ArWei)^2) * sqrt(sum(EKTk.^2 .* ArWei) - sum(EKTk .* ArWei)^2);
    Rho = A / B;
else
    % Unweighted error metrics
    Er = (Diak - EKT);
    SqEr = (Diak - EKT).^2;

    disp('Mean Error');
    ME = mean(Er(~isnan(Er)));

    disp('Mean Absolute Error');
    MAE = mean(abs(Er(~isnan(Er))));

    disp('Root Mean Square Error');
    RMSE = sqrt(mean(SqEr(~isnan(SqEr))));

    disp('Correlation coefficient');
    Rho = corr(EKTk, Diakk);
end

disp('Maximum Absolute Error');
MaxAE = max(abs(Er(~isnan(Er))));

% Compile error metrics into output array
Err = [ME, MAE, RMSE, MaxAE, Rho];
end
