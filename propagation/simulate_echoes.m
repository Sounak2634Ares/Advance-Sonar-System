function rx_data = simulate_echoes(tx_data, environment, targets, beam_config, config)

%% =========================================================
% ECHO SIMULATION MODULE
%
% Purpose:
% Simulates sonar echoes from underwater targets.
%
% Current Features:
% - Target reflections
% - Round-trip delay
% - Echo attenuation
% - Complex noise
%
%% =========================================================

disp('Simulating Echo Returns...');

%% =========================================================
% LOAD PARAMETERS
%% =========================================================

c = config.sound_speed;
fs = config.fs;

tx_signal = tx_data.signal;

signal_length = length(tx_signal);

%% =========================================================
% INITIALIZE RECEIVED SIGNAL
%% =========================================================

rx_signal = zeros(1, signal_length);

%% =========================================================
% TARGET LOOP
%% =========================================================

for k = 1:length(targets)

    % Target position
    x = targets(k).x;
    y = targets(k).y;
    z = targets(k).z;

    % Target reflectivity
    rcs = targets(k).rcs;

    %% -----------------------------------------------------
    % Distance from sonar
    %% -----------------------------------------------------

    distance = sqrt(x^2 + y^2 + z^2);

    %% -----------------------------------------------------
    % Round-trip delay
    %% -----------------------------------------------------

    delay_time = (2 * distance) / c;

    delay_samples = round(delay_time * fs);

    %% -----------------------------------------------------
    % Simple attenuation
    %% -----------------------------------------------------

    attenuation = rcs / (distance^2 + 1);

    %% -----------------------------------------------------
    % Add delayed echo
    %% -----------------------------------------------------

    start_idx = delay_samples + 1;

    if start_idx <= signal_length

        remaining_length = min( ...
            length(tx_signal), ...
            signal_length - delay_samples);

        end_idx = start_idx + remaining_length - 1;

        rx_signal(start_idx:end_idx) = ...
            rx_signal(start_idx:end_idx) + ...
            attenuation .* tx_signal(1:remaining_length);

    end

end

%% =========================================================
% ADD COMPLEX NOISE
%% =========================================================

noise = sqrt(config.noise_power) .* ...
    (randn(size(rx_signal)) + ...
    1i .* randn(size(rx_signal)));

rx_signal = rx_signal + noise;

%% =========================================================
% OUTPUT STRUCTURE
%% =========================================================

rx_data = struct();

rx_data.signal = rx_signal;
rx_data.targets = targets;
rx_data.environment = environment;
rx_data.beam_config = beam_config;

%% =========================================================
% DISPLAY INFORMATION
%% =========================================================

disp('-----------------------------------');
disp('ECHO SIMULATION COMPLETE');
disp('-----------------------------------');

fprintf('Targets Simulated : %d\n', length(targets));
fprintf('Signal Length     : %d samples\n', length(rx_signal));

disp('-----------------------------------');

end