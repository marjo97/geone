function [gssrf] = SSRF_variogram(r, Param)
% SSRF_variogram - Calculates the variance of a stationary spatial random field (SSRF).
%
% This function computes the variance of a stationary spatial random field (SSRF) for a given set
% of distances using specified model parameters. The variance calculation is based on the covariance
% function values at a very small distance (considered as zero lag) and at specified distances.
%
% The SSRF variance is given by:
%   gssrf = varsrf - variosrf + nugget
% where:
%   - varsrf is the covariance at a very small distance (eps, approximated to zero lag),
%   - variosrf is the covariance at the provided distances r,
%   - nugget is an additional term to account for microscale variability or measurement error.
%
% INPUTS:
%   r     - A vector of distances at which to calculate the SSRF variance.
%   Param - A vector of model parameters for the SSRF covariance function [s^2, b, c, nugget], where:
%           s^2: Variance component.
%           b  : Scale parameter, controlling how quickly covariance decays with distance.
%           c  : Shape parameter (not utilized directly in this function).
%           nugget: Nugget effect, accounting for discontinuities or micro-scale variations.
%
% OUTPUT:
%   gssrf - The computed SSRF variance values corresponding to each distance in r.
%
% EXAMPLE USAGE:
%   Param = [1.0, 2.0, 1.195, 0.1]; % Example SSRF model parameters
%   r = linspace(0, 5, 100);        % Define a range of distances
%   gssrf = SSRF_variogram(r, Param); % Calculate SSRF variance
%   plot(r, gssrf);                 % Plot the variance as a function of distance

% Calculate the covariance at a very small distance (close to zero lag) using cov_ssrf_exp
[~, varsrf] = SSRF_cov(3, eps*1000, Param, 0);

% Calculate the covariance at the specified distances r
[~, variosrf] = SSRF_cov(3, r, Param, 0);

% Compute the SSRF variance by subtracting covariance values and adding the nugget effect
gssrf = varsrf - variosrf; % Difference between zero lag and distance r covariance
gssrf = gssrf + Param(4);  % Add nugget effect to account for measurement error or microscale variation

end

