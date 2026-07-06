function ui_state = render_visualization(beam_energy, theta_vec, range_vec, detectionsXYZ, environment, ui_state)
% RENDER_VISUALIZATION High-performance industrial sonar waterfall UI renderer.
%
% File        : render_visualization.m
% Purpose     : Phase 3c Graphic Matrix Visualization Shading Engine

arguments
    beam_energy {mustBeNumeric, mustBeNonempty}
    theta_vec {mustBeNumeric}
    range_vec = []
    detectionsXYZ = []
    environment struct = struct()
    ui_state = struct()
end

% Configuration fallback initialization
if ~isfield(environment, 'dynamicRange_dB'), environment.dynamicRange_dB = 60; end
if ~isfield(environment, 'gainDisplay'),       environment.gainDisplay = 1.0; end

[NumSamples, NumBeams] = size(beam_energy);

% Synthesize arrays configurations if tracking vectors are left empty
if isempty(range_vec), range_vec = (0:NumSamples-1).'; end
if isempty(theta_vec), theta_vec = linspace(-pi/3, pi/3, NumBeams); end

% Clean up orientations
theta_vec = theta_vec(:).'; 
range_vec = range_vec(:);   

% Establish axes hooks bindings securely
if ~isfield(ui_state, 'axesWaterfall') || isempty(ui_state.axesWaterfall)
    figure('Name', 'Sonar Visualization Port', 'Color', [0.05 0.06 0.08]);
    ax1 = subplot(1, 2, 1); ax2 = subplot(1, 2, 2);
    ui_state.axesWaterfall = ax1;
    ui_state.axes3D = ax2;
end

axPolar = ui_state.axesWaterfall;
ax3d = ui_state.axes3D;

% Convert raw power to high-dynamic-range log scale decibel presentation levels
dB_map = 10 * log10(beam_energy + eps);
dB_map = dB_map * environment.gainDisplay;

% Calculate bounding display limits clamping values
max_val = max(dB_map(:));
min_clamp = max_val - environment.dynamicRange_dB;

% =========================================================================
% RENDER POLAR WATERFALL (SURF-TEXTURE SHADERS GRAPHICS PIPELINE)
% =========================================================================
% Construct polar mesh grid matrices shapes
[ThetaMesh, RangeMesh] = meshgrid(theta_vec, range_vec);
X_mesh = RangeMesh .* sin(ThetaMesh);
Y_mesh = RangeMesh .* cos(ThetaMesh);
Z_mesh = zeros(size(X_mesh));

if ~isfield(ui_state, 'h_surf') || ~ishandle(ui_state.h_surf)
    cla(axPolar);
    ui_state.h_surf = surf(axPolar, X_mesh, Y_mesh, Z_mesh, dB_map, ...
        'EdgeColor', 'none', 'FaceColor', 'texturemap');
    % 'cyan' is not a built-in MATLAB colormap name (valid built-ins are
    % parula, turbo, hot, cool, etc.) -- this line previously threw
    % "Unrecognized function or variable 'cyan'" and aborted rendering.
    % Built a genuine black -> cyan gradient to match the intended
    % waterfall aesthetic instead of guessing at a named colormap.
    cyan_cmap = [linspace(0, 0.05, 256).', linspace(0, 0.85, 256).', linspace(0.05, 1, 256).'];
    colormap(axPolar, cyan_cmap);
    view(axPolar, 2);
    axis(axPolar, 'equal', 'tight');
    grid(axPolar, 'on');
    axPolar.Color = [0 0 0];
    axPolar.XColor = [0.5 0.7 0.9];
    axPolar.YColor = [0.5 0.7 0.9];
    title(axPolar, 'Fast-Time Polar Beam Intensity Map', 'Color', 'w');
else
    set(ui_state.h_surf, 'XData', X_mesh, 'YData', Y_mesh, 'CData', dB_map);
end
clim(axPolar, [min_clamp, max_val]);

% =========================================================================
% RENDER 3D TARGET OVERLAY INTERFACE DISPLAY
% =========================================================================
if ~isempty(detectionsXYZ)
    X_det = detectionsXYZ(:, 1);
    Y_det = detectionsXYZ(:, 2);
    Z_det = detectionsXYZ(:, 3);
else
    X_det = []; Y_det = []; Z_det = [];
end

if ~isfield(ui_state, 'h_scatter') || ~ishandle(ui_state.h_scatter)
    cla(ax3d);
    ui_state.h_scatter = scatter3(ax3d, X_det, Y_det, Z_det, 50, 'filled', ...
        'MarkerFaceColor', [1 0.2 0.2], 'MarkerEdgeColor', [1 0.8 0.8]);
    grid(ax3d, 'on');
    view(ax3d, 3);
    ax3d.Color = [0 0 0];
    ax3d.XColor = [0.5 0.7 0.9];
    ax3d.YColor = [0.5 0.7 0.9];
    ax3d.ZColor = [0.5 0.7 0.9];
    xlabel(ax3d, 'Cross-Range X (m)'); ylabel(ax3d, 'Boresight Y (m)'); zlabel(ax3d, 'Depth Z (m)');
    title(ax3d, 'Reconstructed Target Operations Space', 'Color', 'w');
else
    set(ui_state.h_scatter, 'XData', X_det, 'YData', Y_det, 'ZData', Z_det);
end

drawnow limitrate;

end