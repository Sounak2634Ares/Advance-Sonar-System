function tx_signal = generate_chirp(config)

%% =========================================================
% SONAR CHIRP GENERATION MODULE
%
% Purpose:
% Generates the transmitted sonar waveform.
%
% Features:
% - Linear Frequency Modulated (LFM) Chirp
% - Adjustable bandwidth
% - Adjustable pulse duration
% - Windowing support
% - Future GPU compatibility
%% =========================================================

disp('Generating Sonar Chirp Signal...');

%% =========================================================
% LOAD PARAMETERS FROM CONFIG
%% =========================================================

fc = config.fc;

bandwidth = config.bandwidth;

pulse_duration = config.pulse_duration;

fs = config.fs;

%% =========================================================
% TIME VECTOR
%% =========================================================

t = 0 : 1/fs : pulse_duration;

%% =========================================================
% CHIRP RATE
%% =========================================================

k = bandwidth / pulse_duration;

%% =========================================================
% GENERATE LFM CHIRP
%% =========================================================

tx_signal = exp( ...
            1j * 2 * pi * ...
            ( ...
            fc * t ...
            + 0.5 * k * t.^2 ...
            ) ...
            );

%% =========================================================
% APPLY WINDOW FUNCTION
%% =========================================================

window = hamming(length(tx_signal))';

tx_signal = tx_signal .* window;

%% =========================================================
% NORMALIZE SIGNAL
%% =========================================================

tx_signal = tx_signal ./ max(abs(tx_signal));

%% =========================================================
% DISPLAY SIGNAL INFORMATION
%% =========================================================

disp('-----------------------------------');
disp('CHIRP SIGNAL GENERATED');
disp('-----------------------------------');

fprintf('Center Frequency : %.2f Hz\n', fc);

fprintf('Bandwidth        : %.2f Hz\n', bandwidth);

fprintf('Pulse Duration   : %.5f s\n', pulse_duration);

fprintf('Sampling Rate    : %.2f Hz\n', fs);

fprintf('Signal Length    : %d samples\n', length(tx_signal));

disp('-----------------------------------');

end