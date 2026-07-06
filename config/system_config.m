function config = system_config()

%% =========================================================
% SYSTEM CONFIGURATION FILE
%
% Purpose:
% Stores all global sonar parameters.
%
% Every module will use this configuration structure.
%% =========================================================

disp('Loading Sonar System Configuration...');

%% =========================================================
% GENERAL SYSTEM PARAMETERS
%% =========================================================

config.system_name = 'Advanced Multi-Beam Sonar System';

config.version = '1.0';

config.author = 'Sounak Jana';

%% =========================================================
% SONAR PARAMETERS
%% =========================================================

%% FIX:
% Every downstream module (signal_processing.m, beamforming.m,
% target_detection.m, generate_sonar_echoes.m, reconstruct_3d.m) reads
% config.Fs, config.fc, config.BW, config.T, config.N, config.d, config.c,
% config.NumBeams, config.beamformer_type, config.pfa (capitalized, short
% form). This file previously only produced lowercase/long-form fields
% (fs, bandwidth, sound_speed, num_beams, array_elements, element_spacing)
% that no other module ever reads, so calling system_config() and passing
% its output into the pipeline would silently fall back to hard-coded
% defaults everywhere. Canonical pipeline-facing fields are added below;
% the original descriptive fields are kept as aliases for compatibility.

config.fc = 50000;              % Center frequency (Hz)

config.bandwidth = 60000;       % Legacy alias, kept for compatibility
config.BW = config.bandwidth;   % Bandwidth (Hz) -- canonical name

config.pulse_duration = 0.01;
config.T = config.pulse_duration; % Pulse duration (s) -- canonical name

config.fs = 1000000;
config.Fs = config.fs;          % Sampling rate (Hz) -- canonical name

config.sound_speed = 1500;
config.c = config.sound_speed;  % Speed of sound (m/s) -- canonical name

config.SNR_dB = 20;              % Operational SNR target (not present before)

config.beamformer_type = 'MVDR'; % Choice of spatial processor: 'DAS' or 'MVDR'

config.pfa = 1e-4;               % CFAR constant false-alarm-rate bound

%% =========================================================
% MULTI-BEAM CONFIGURATION
%% =========================================================

config.num_beams = 64;
config.NumBeams = config.num_beams; % canonical name used by beamforming.m

config.max_beam_angle = 60;

config.array_elements = 32;
config.N = config.array_elements;   % canonical name used by beamforming.m

config.element_spacing = 0.015;
config.d = config.element_spacing;  % canonical name used by beamforming.m

%% =========================================================
% ENVIRONMENT SETTINGS
%% =========================================================

config.max_range = 100;

config.max_depth = 50;

config.water_temperature = 25;

config.salinity = 35;

%% =========================================================
% NOISE SETTINGS
%% =========================================================

config.noise_power = 0.01;

config.reverberation_level = 0.02;

%% =========================================================
% GPU SETTINGS
%% =========================================================

config.use_gpu = true;

%% =========================================================
% VISUALIZATION SETTINGS
%% =========================================================

config.enable_3d = true;

config.enable_waterfall = true;

config.enable_real_time = true;

%% =========================================================
% DATA STORAGE SETTINGS
%% =========================================================

config.save_results = true;

config.output_folder = 'data';

%% =========================================================
% CONFIGURATION COMPLETE
%% =========================================================

disp('System Configuration Loaded Successfully.');

end