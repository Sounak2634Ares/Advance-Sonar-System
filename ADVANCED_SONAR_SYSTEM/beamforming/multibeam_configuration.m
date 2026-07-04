function beam_config = multibeam_configuration(config)

%% =========================================================
% MULTI-BEAM SONAR CONFIGURATION MODULE
%
% Purpose:
% Configures sonar array geometry and beam steering.
%
% Features:
% - Uniform Linear Array (ULA)
% - Beam steering angles
% - Element positions
% - Future beamforming support
%
%% =========================================================

disp('Configuring Multi-Beam Sonar Array...');

%% =========================================================
% LOAD PARAMETERS
%% =========================================================

num_beams = config.num_beams;

max_beam_angle = config.max_beam_angle;

array_elements = config.array_elements;

element_spacing = config.element_spacing;

%% =========================================================
% BEAM ANGLES
%% =========================================================

beam_angles = linspace( ...
                -max_beam_angle, ...
                 max_beam_angle, ...
                 num_beams);

%% =========================================================
% ARRAY GEOMETRY
%% =========================================================

element_positions = ...
    (0:array_elements-1) * element_spacing;

element_positions = ...
    element_positions - mean(element_positions);

%% =========================================================
% STEERING MATRIX
%% =========================================================

steering_matrix = zeros( ...
                    array_elements, ...
                    num_beams);

for beam = 1:num_beams

    theta = deg2rad(beam_angles(beam));

    steering_matrix(:,beam) = ...
        exp(-1j * 2*pi * ...
        element_positions' * sin(theta));

end

%% =========================================================
% STORE CONFIGURATION
%% =========================================================

beam_config = struct();

beam_config.num_beams = num_beams;

beam_config.max_beam_angle = max_beam_angle;

beam_config.array_elements = array_elements;

beam_config.element_spacing = element_spacing;

beam_config.beam_angles = beam_angles;

beam_config.element_positions = element_positions;

beam_config.steering_matrix = steering_matrix;

%% =========================================================
% DISPLAY INFORMATION
%% =========================================================

disp('-----------------------------------');
disp('MULTI-BEAM CONFIGURATION COMPLETE');
disp('-----------------------------------');

fprintf('Number of Beams      : %d\n', num_beams);

fprintf('Array Elements       : %d\n', array_elements);

fprintf('Element Spacing      : %.4f m\n', ...
        element_spacing);

fprintf('Beam Coverage        : %.1f° to %.1f°\n', ...
        min(beam_angles), ...
        max(beam_angles));

disp('-----------------------------------');

end