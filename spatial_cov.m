function [C] = spatial_cov(Param, r)
% spatial_cov - Computes spatial covariance using a specific model.
%
% This function calculates the spatial covariance for a two-dimensional field based
% on a specific covariance model outlined in Hristopoulos' paper (Equation 46a).
% The model is suitable for isotropic conditions and includes a nugget effect
% to handle measurement errors or micro-scale variations.
%
% The covariance formula for d = 2 (two-dimensional case) is:
%   C = s^2 * 0.1592 * b / (r^2 + b^2)^(3/2)
% where:
%   - s^2 is the variance component (V),
%   - b is a scale parameter,
%   - r is the distance,
%   - nugget is the nugget effect to model discontinuities at the origin (small-scale variability).
%
% INPUTS:
%   Param - A vector containing the model parameters [s^2, b, nugget]
%           s^2    : Variance component, representing the sill of the variogram.
%           b      : Scale parameter, controlling the range over which the covariance diminishes.
%           nugget : Nugget effect, accounting for measurement error or microscale variation.
%   r     - A vector of distances at which to calculate the covariance.
%
% OUTPUT:
%   C - The computed covariance values for each distance in r.
%
% EXAMPLE USAGE:
%   Param = [1.0, 2.0, 0.1]; % Example parameters
%   r = linspace(0, 5, 100); % Distance range from 0 to 5
%   C = spatial_cov(Param, r); % Compute covariance
%   plot(r, C); % Plot the covariance function
%
% NOTES:
%   - This function is designed for 2D spatial analysis.
%   - The nugget effect is added to handle cases where the distance r is effectively zero (r <= eps).
%   - The constant 0.1592 is derived from the model specification for scaling the covariance function.


% Extract parameters from input
V = Param(1); % V is the variance, s^2
b = Param(2); % Scale parameter
nug = Param(3); % Nugget effect

% Compute the covariance using the specified model formula
C = V * 0.1592 * b ./ ((r.^2 + b^2).^(3/2)) + (nug * (r <= eps));

end

