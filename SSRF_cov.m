function [h, covar] = SSRF_cov(dd, lag, Param, flag)
% SSRF_cov - Calculates the explicit Spartan Spectral Random Field (SSRF) covariance function.
%
% This function computes the explicit SSRF covariance function, which is useful in
% geostatistics for modeling spatial dependencies. The model includes a nugget effect
% and supports various dimensionalities (1D, 2D, 3D). The covariance is calculated
% at the limit of an infinite spectral cutoff, allowing for flexible modeling of spatial
% data with different ranges and smoothness.
%
% INPUTS:
%   dd    - Dimensionality of the field (1, 2, 3, etc.). Determines the spatial domain.
%   lag   - Vector of lag distances over which to compute the covariance.
%   Param - Vector of Spartan parameters [eta0, eta1, xi, nugget]:
%           * eta0: Amplitude parameter of the covariance function.
%           * eta1: Shape parameter affecting the decay rate.
%           * xi:   Range parameter influencing the spatial correlation length.
%           * nugget: Nugget effect representing microscale variability (optional).
%   flag  - Plot flag (0 = suppress plots, 1 = show plots).
%
% OUTPUTS:
%   h     - Normalized lag vector (lag distances divided by range parameter xi).
%   covar - Covariance vector calculated for each lag distance.
%
% EXAMPLE USAGE:
%   [h, covar] = SSRF_cov(2, 0:1:10, [1, 2, 3, 0.5], 1);
%
% REFERENCES:
%   D. T. Hristopulos, Spartan Gibbs random field models for geostatistical
%   applications, SIAM Journal of Scientific Computing, 24 (6), pp. 2125-2162 (2003).
%
% NOTE:
%   This function has been tested and is compatible with both Octave and MATLAB environments.
%   It provides different covariance models based on the dimensionality and parameters.

% Extract parameters from input
eta0 = Param(1);
eta1 = Param(2);
xi = Param(3);
nparam = length(Param);

% Calculate the normalized lag vector
h = lag / xi;

% Replace zero lag with a small value to avoid division by zero (De L'Hopital limit)
h(lag == 0) = eps;

% Initialize covariance function based on dimensionality (dd)
switch dd
    case 1 % 1D Covariance
        if abs(eta1) < 2
            beta1 = 0.5 * sqrt(2 - eta1);
            beta2 = 0.5 * sqrt(2 + eta1);
            Gr = (exp(-h * beta2)) .* ((cos(h * beta1)) / beta2 + (sin(h * beta1)) / beta1) / 4;
        elseif eta1 == 2
            Gr = ((1 + h) .* exp(-h)) / 4;
        elseif eta1 > 2
            delta = sqrt(eta1^2 - 4);
            omega2 = sqrt((eta1 + delta) / 2);
            omega1 = sqrt((eta1 - delta) / 2);
            Gr = (exp(-h * omega1) / omega1 - exp(-h * omega2) / omega2) / (2 * delta);
        else
            disp('Permissibility conditions fall');
            Gr = nan(size(h));
        end
    case 2 % 2D Covariance
        if eta1 == 2
            Gr = h .* besselk(-1, h) / (4 * pi);
            vv = 1 / (4 * pi);
        elseif eta1 > 2
            Dn = sqrt(eta1^2 - 4);
            zp = sqrt(-(-eta1 + Dn) / 2);
            zm = sqrt(-(-eta1 - Dn) / 2);
            Gr = (besselk(0, h * zp) - besselk(0, h * zm)) / (2 * pi * Dn);
            vv = log((eta1 + Dn) / (eta1 - Dn)) / (4 * pi * Dn);
        elseif abs(eta1) < 2
            Dn = sqrt(4 - eta1^2);
            zp = sqrt(-(-eta1 + sqrt(-1) * Dn) / 2);
            zm = sqrt(-(-eta1 - sqrt(-1) * Dn) / 2);
            Gr = imag(besselk(0, h * zp)) / (pi * Dn);
            vv = (pi / 2 - atan(eta1 / Dn)) / (2 * pi * Dn);
        else
            disp('Permissibility conditions fall');
            Gr = nan(size(h));
        end
    case 3 % 3D Covariance
        if abs(eta1) < 2
            delta = sqrt(4 - eta1^2);
            beta1 = 0.5 * sqrt(2 - eta1);
            beta2 = 0.5 * sqrt(2 + eta1);
            Gr = (exp(-h * beta2)) .* (sin(h * beta1) ./ (2 * pi * h * delta));
        elseif eta1 == 2
            Gr = exp(-h) / (8 * pi);
        elseif eta1 > 2
            delta = sqrt(eta1^2 - 4);
            omega2 = sqrt((eta1 + delta) / 2);
            omega1 = sqrt((eta1 - delta) / 2);
            Gr = (exp(-h * omega1) - exp(-h * omega2)) ./ (4 * pi * h * delta);
        else
            disp('Permissibility conditions fall'); % Handle invalid eta1 values
            Gr = nan(size(h));
        end
    otherwise
        disp('Not defined for d > 3'); % Handle unsupported dimensionalities
        Gr = nan(size(h));
end

% Calculate the covariance, incorporating the nugget effect if specified
if eta1 <= -2
    covar = nan(size(h)); % Invalid eta1 results in NaN covariance
else
    if nparam == 4
        C0 = Param(4); % Nugget effect
        covar = eta0 * Gr + C0;
    else
        covar = eta0 * Gr; % Without nugget effect
    end
end

% Plotting the covariance if the flag is set
if nargin > 3 && flag > 0
    figure;
    plot(h, covar, 'b-'); % Plot covariance curve
    if length(h) < 20
        hold on;
        plot(h, covar, 'r.'); % Highlight individual points if few
        hold off;
    end
    xlabel('h');
    ylabel('SRRF Covariance');
    title(['d=', num2str(dd)]);
    axis tight;
end

end

