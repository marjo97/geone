function Apo = kriging_after_bt(Xk, Yk, VAL, SF)
% Krigrids_after_bt - Displays Kriging maps and computes volume estimates.
%
% This function creates two types of plots for Kriging results: a surface plot
% of the predicted values and an optional surface plot of the Kriging standard
% deviation. It also calculates the volume of the deposit in the same units as
% the input coordinates.
%
% INPUTS:
%   Xk, Yk - Vectors of X and Y coordinates for the grid points
%   VAL    - Matrix of estimated values at the prediction points (fluctuation or random field values)
%   SF     - Matrix of the Kriging standard deviation at the prediction points (optional)
%
% OUTPUT:
%   Apo    - The estimated volume of the deposit in the same units as Xk and Yk
%
% NOTE: The function retains NaN values for predictions where not enough neighbors are available.
%
% EXAMPLE USAGE:
%   Xk = linspace(0, 10, 50);
%   Yk = linspace(0, 10, 50);
%   [Xg, Yg] = meshgrid(Xk, Yk);
%   VAL = sin(Xg) + cos(Yg);
%   SF = rand(size(Xg)) * 0.1;
%   Apo = Krigrids_after_bt(Xk, Yk, VAL, SF);

%% Figure of Predictions
% Create a surface plot for the predicted values
figure;
surf(Xk, Yk, VAL, 'EdgeColor', 'none');
shading flat;
colorbar;
colormap('jet');
xlabel('West-East', 'FontSize', 25, 'FontWeight', 'bold');
ylabel('South-North', 'FontSize', 25, 'FontWeight', 'bold');
set(gca, 'FontSize', 30, 'LineWidth', 1.5);
grid on;
axis equal tight;
view(2);

%% Figure of Kriging error
% If the Kriging standard deviation matrix SF is provided, create a plot for it
if nargin > 3
    figure;
    surf(Xk, Yk, SF, 'EdgeColor', 'none');
    shading flat;
    colorbar;
    colormap('bone');
    title('Kriging error standard deviation', 'FontSize', 20, 'FontWeight', 'bold');
    xlabel('West-East', 'FontSize', 18, 'FontWeight', 'bold');
    ylabel('South-North', 'FontSize', 18, 'FontWeight', 'bold');
    set(gca, 'FontSize', 18, 'LineWidth', 1.5);
    grid on;
    axis equal tight;
    view(2);
end

%% Sum of the Prediction
% Calculate the volume estimate by summing up the predicted values
% Calculate the area of each grid cell based on the spacing in X and Y directions
CellX = Xk(2) - Xk(1);  % Spacing in the X direction
CellY = Yk(2) - Yk(1);  % Spacing in the Y direction
ArCell = CellX * CellY; % Area of each cell

% Calculate the estimated volume (Apo) by summing up all predicted values, ignoring NaNs
Apo = ArCell * nansum(nansum(VAL)); % Sum up the values, ignoring NaNs
end

