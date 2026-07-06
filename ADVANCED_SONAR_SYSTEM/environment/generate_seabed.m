function environment = generate_seabed(config)

%% =========================================================
% SEABED GENERATION MODULE
%
% Purpose:
% Creates underwater terrain for sonar simulation.
%
% Features:
% - Random seabed generation
% - Terrain roughness
% - Depth map generation
% - Future texture support
%% =========================================================

disp('Generating Seabed Environment...');

%% =========================================================
% LOAD PARAMETERS
%% =========================================================

max_range = config.max_range;

max_depth = config.max_depth;

%% =========================================================
% CREATE GRID
%% =========================================================

x = linspace(0, max_range, 500);

y = linspace(-50, 50, 300);

[X, Y] = meshgrid(x, y);

%% =========================================================
% BASE SEABED PROFILE
%% =========================================================

Z = -20 ...
    - 5 * sin(0.05 * X) ...
    - 3 * cos(0.08 * Y);

%% =========================================================
% ADD RANDOM ROUGHNESS
%% =========================================================

roughness = 0.8 * randn(size(Z));

Z = Z + roughness;

%% =========================================================
% LIMIT MAXIMUM DEPTH
%% =========================================================

Z(Z < -max_depth) = -max_depth;

%% =========================================================
% STORE ENVIRONMENT DATA
%% =========================================================

environment.X = X;

environment.Y = Y;

environment.Z = Z;

environment.roughness = roughness;

%% =========================================================
% DISPLAY INFORMATION
%% =========================================================

disp('-----------------------------------');
disp('SEABED GENERATED');
disp('-----------------------------------');

fprintf('Environment Range : %.2f meters\n', max_range);

fprintf('Maximum Depth     : %.2f meters\n', max_depth);

fprintf('Grid Size         : %d x %d\n', ...
        size(X,1), size(X,2));

disp('-----------------------------------');

end