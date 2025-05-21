function combinedEigenfunctionPlot(X, Z, Xi, Psi, B)
% combinedEigenfunctionPlot creates a publication-quality figure that contains
% both the full 2D eigenfunction plots and a selected x-slice (midpoint) for each
% eigenfunction. The top row shows the full-field plots using pcolor and the bottom
% row shows the eigenfunction values along a fixed x location.
%
%   combinedEigenfunctionPlot(X, Z, Xi, Psi, B) uses:
%       X, Z  - 2D arrays defining the spatial grid (as from meshgrid).
%       Xi    - 2D eigenfunction field (first field).
%       Psi   - 2D eigenfunction field (second field).
%       B     - 2D eigenfunction field (third field).
%
%   The top row displays the fields over (x,z) and the bottom row displays
%   the eigenfunction profiles along a specific x value (midpoint in x).
%
%   Example:
%       [X, Z] = meshgrid(linspace(0,10,100), linspace(0,5,50));
%       Xi = sin(X) + cos(Z);
%       Psi = X.^2 - Z.^2;
%       B = exp(-((X-5).^2+(Z-2.5).^2));
%       combinedEigenfunctionPlot(X, Z, Xi, Psi, B)

%% Aesthetic settings
fs = 16;    % Font size for labels and titles
lw = 3;     % Line width for plots and axes
cmap = turbo;  % Use the "turbo" colormap

% Set default interpreter for LaTeX on axes tick labels, legends, etc.
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');

%% Create a full-screen figure with a tiled layout.
fig = figure('WindowState', 'maximized');
% Create a tiled layout with 2 rows and 3 columns.
tiledlayout(2,3, 'Padding', 'compact', 'TileSpacing','compact');

%% ------------------------------
%% Top Row: Full Eigenfunction Plots using pcolor
%% ------------------------------

% --- Tile 1: Field Xi ---
nexttile;
pcolor(X, Z, Xi);
shading interp;     % smooth color transitions
axis tight;
daspect([1 1 1]);   % set data aspect ratio 1:1
colormap(cmap);
c = colorbar;
c.Label.Interpreter = 'latex';
c.Label.FontSize = fs;
c.Label.String = '$\xi$';
xlabel('$x/Fr$', 'Interpreter', 'latex', 'FontSize', fs);
ylabel('$z$', 'Interpreter', 'latex', 'FontSize', fs);
title('Eigenfunction $\xi$', 'Interpreter', 'latex', 'FontSize', fs);
box on;
set(gca, 'LineWidth', lw, 'FontSize', fs, 'BoxStyle', 'full');

% --- Tile 2: Field Psi ---
nexttile;
pcolor(X, Z, Psi);
shading interp;
axis tight;
daspect([1 1 1]);
colormap(cmap);
c = colorbar;
c.Label.Interpreter = 'latex';
c.Label.FontSize = fs;
c.Label.String = '$\psi$';
xlabel('$x/Fr$', 'Interpreter', 'latex', 'FontSize', fs);
ylabel('$z$', 'Interpreter', 'latex', 'FontSize', fs);
title('Eigenfunction $\psi$', 'Interpreter', 'latex', 'FontSize', fs);
box on;
set(gca, 'LineWidth', lw, 'FontSize', fs, 'BoxStyle', 'full');

% --- Tile 3: Field B ---
nexttile;
pcolor(X, Z, B);
shading interp;
axis tight;
daspect([1 1 1]);
colormap(cmap);
c = colorbar;
c.Label.Interpreter = 'latex';
c.Label.FontSize = fs;
c.Label.String = '$b$';
xlabel('$x/Fr$', 'Interpreter', 'latex', 'FontSize', fs);
ylabel('$z$', 'Interpreter', 'latex', 'FontSize', fs);
title('Eigenfunction $b$', 'Interpreter', 'latex', 'FontSize', fs);
box on;
set(gca, 'LineWidth', lw, 'FontSize', fs, 'BoxStyle', 'full');

%% ------------------------------
%% Bottom Row: x-Slice Plots (midpoint of x)
%% ------------------------------

% Determine the mid-point index along the x-direction.
Nx = size(X,2);
midIdx = floor(0.5 * Nx) + 1;

% Extract the x-slice from each field (all rows, fixed midIdx column).
xi_slice  = Xi(:, midIdx);
psi_slice = Psi(:, midIdx);
b_slice   = B(:, midIdx);

% Since meshgrid makes every column of Z identical, extract the z coordinate from the first column.
z_vec = Z(:,1);

% --- Tile 4: Slice for Xi ---
nexttile;
plot(xi_slice, z_vec, 'LineWidth', lw);
grid on;
box on;
xlabel('Slice $\xi$', 'Interpreter', 'latex', 'FontSize', fs);
ylabel('$z$', 'Interpreter', 'latex', 'FontSize', fs);
title('Eigenfunction Slice $\xi$', 'Interpreter', 'latex', 'FontSize', fs);
set(gca, 'LineWidth', lw, 'FontSize', fs, 'BoxStyle', 'full');
axis tight;

% --- Tile 5: Slice for Psi ---
nexttile;
plot(psi_slice, z_vec, 'LineWidth', lw);
grid on;
box on;
xlabel('Slice $\psi$', 'Interpreter', 'latex', 'FontSize', fs);
ylabel('$z$', 'Interpreter', 'latex', 'FontSize', fs);
title('Eigenfunction Slice $\psi$', 'Interpreter', 'latex', 'FontSize', fs);
set(gca, 'LineWidth', lw, 'FontSize', fs, 'BoxStyle', 'full');
axis tight;

% --- Tile 6: Slice for B ---
nexttile;
plot(b_slice, z_vec, 'LineWidth', lw);
grid on;
box on;
xlabel('Slice $b$', 'Interpreter', 'latex', 'FontSize', fs);
ylabel('$z$', 'Interpreter', 'latex', 'FontSize', fs);
title('Eigenfunction Slice $b$', 'Interpreter', 'latex', 'FontSize', fs);
set(gca, 'LineWidth', lw, 'FontSize', fs, 'BoxStyle', 'full');
axis tight;

end
