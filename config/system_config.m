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

config.fc = 50000;

config.bandwidth = 60000;

config.pulse_duration = 0.01;

config.fs = 1000000;

config.sound_speed = 1500;

%% =========================================================
% MULTI-BEAM CONFIGURATION
%% =========================================================

config.num_beams = 64;

config.max_beam_angle = 60;

config.array_elements = 32;

config.element_spacing = 0.015;

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