function tx_data = transmit_signal( ...
                        tx_signal, ...
                        environment, ...
                        beam_config, ...
                        config)

%% =========================================================
% SIGNAL TRANSMISSION MODULE
%
% Purpose:
% Simulates transmission of sonar pulse through water.
%
% Features:
% - Multi-beam transmission
% - Beam angle assignment
% - Range propagation setup
% - Future transmission loss support
%
%% =========================================================

disp('Transmitting Sonar Signal...');

%% =========================================================
% LOAD PARAMETERS
%% =========================================================

sound_speed = config.sound_speed;

max_range = config.max_range;

num_beams = beam_config.num_beams;

beam_angles = beam_config.beam_angles;

%% =========================================================
% CREATE TRANSMISSION STRUCTURE
%% =========================================================

tx_data = struct();

tx_data.signal = tx_signal;

tx_data.sound_speed = sound_speed;

tx_data.max_range = max_range;

tx_data.num_beams = num_beams;

tx_data.beam_angles = beam_angles;

tx_data.environment = environment;

%% =========================================================
% CALCULATE RANGE AXIS
%% =========================================================

range_axis = linspace( ...
                0, ...
                max_range, ...
                length(tx_signal));

tx_data.range_axis = range_axis;

%% =========================================================
% CALCULATE BEAM DIRECTIONS
%% =========================================================

beam_vectors = zeros(num_beams,2);

for k = 1:num_beams

    theta = deg2rad(beam_angles(k));

    beam_vectors(k,1) = cos(theta);
    beam_vectors(k,2) = sin(theta);

end

tx_data.beam_vectors = beam_vectors;

%% =========================================================
% DISPLAY INFORMATION
%% =========================================================

disp('-----------------------------------');
disp('TRANSMISSION CONFIGURED');
disp('-----------------------------------');

fprintf('Number of Beams : %d\n', num_beams);

fprintf('Maximum Range   : %.2f m\n', max_range);

fprintf('Sound Speed     : %.2f m/s\n', sound_speed);

fprintf('Signal Length   : %d samples\n', ...
        length(tx_signal));

disp('-----------------------------------');

end