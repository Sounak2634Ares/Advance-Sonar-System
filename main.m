%% =========================================================
% ADVANCED MULTI-BEAM SONAR SYSTEM
% MASTER CONTROL FILE
%
% File        : main.m
% Purpose     : Central execution controller
% Author      : Sounak Jana
% MATLAB      : R2024b
%
% Description :
% This file controls the entire sonar framework.
% No processing logic should be written here.
% Only:
%   - Initialization
%   - Module calling
%   - Execution control
%   - Error handling
%   - System flow management
%
% All major processing must remain inside modules.
%% =========================================================

clc;
clear;
close all;

format compact;

%% =========================================================
% PROJECT INITIALIZATION
%% =========================================================

disp('=================================================');
disp(' ADVANCED MULTI-BEAM SONAR SYSTEM INITIALIZING ');
disp('=================================================');

%% =========================================================
% ADD PROJECT PATHS
%% =========================================================

disp('Loading Project Paths...');

addpath(genpath('config'));
addpath(genpath('signal_generation'));
addpath(genpath('environment'));
addpath(genpath('targets'));
addpath(genpath('beamforming'));
addpath(genpath('propagation'));
addpath(genpath('processing'));
addpath(genpath('reconstruction'));
addpath(genpath('visualization'));
addpath(genpath('gpu'));
addpath(genpath('ui'));
addpath(genpath('hardware'));
addpath(genpath('data'));
addpath(genpath('utils'));

disp('Project Paths Loaded.');

%% =========================================================
% START EXECUTION TIMER
%% =========================================================

total_runtime = tic;

%% =========================================================
% SYSTEM CONFIGURATION
%% =========================================================

disp('Loading System Configuration...');

config = system_config();

disp('System Configuration Loaded.');

%% =========================================================
% GPU INITIALIZATION
%% =========================================================

disp('Initializing GPU...');

gpu_status = initialize_gpu(config);

disp('GPU Initialization Complete.');

%% =========================================================
% SONAR SIGNAL GENERATION
%% =========================================================

disp('Generating Sonar Waveform...');

tx_signal = generate_chirp(config);

disp('Waveform Generation Complete.');

%% =========================================================
% UNDERWATER ENVIRONMENT GENERATION
%% =========================================================

disp('Generating Underwater Environment...');

environment = generate_seabed(config);

disp('Environment Generation Complete.');

%% =========================================================
% TARGET GENERATION
%% =========================================================

disp('Generating Targets...');

targets = generate_targets(config);

disp('Target Generation Complete.');

%% =========================================================
% MULTI-BEAM ARRAY CONFIGURATION
%% =========================================================

disp('Configuring Sonar Array...');

beam_config = multibeam_configuration(config);

disp('Beamforming Configuration Complete.');

%% =========================================================
% SIGNAL TRANSMISSION
%% =========================================================

disp('Transmitting Sonar Signal...');

tx_data = transmit_signal( ...
            tx_signal, ...
            environment, ...
            beam_config, ...
            config);

disp('Signal Transmission Complete.');

%% =========================================================
% ECHO SIMULATION
%% =========================================================

disp('Simulating Echo Returns...');

rx_data = simulate_echoes( ...
            tx_data, ...
            environment, ...
            targets, ...
            beam_config, ...
            config);

disp('Echo Simulation Complete.');

%% =========================================================
% SIGNAL PROCESSING
%% =========================================================

disp('Processing Sonar Data...');

% Execute optimized fast-time pulse compression and TVG normalization
[processed_data, proc_meta] = signal_processing( ...
                                rx_data, ...
                                tx_signal, ...
                                config);

% --- DATA GUARD RAILS & INTEGRITY CHECKS ---
assert(isequal(size(processed_data), [config.NumSamples config.N]), ...
    'Dimension Error: signal_processing output array dimensions do not match [NumSamples x N].');

assert(isfloat(processed_data), ...
    'Data Type Error: processed_data must be floating-point values.');

assert(all(isfinite(processed_data(:))), ...
    'Numerical Instability Detected: Matched Filtering or TVG operation generated NaN/Inf values.');

disp('Signal Processing Complete.');

%% =========================================================
% TARGET DETECTION
%% =========================================================

disp('Running Target Detection...');

detections = target_detection( ...
                processed_data, ...
                config);

disp('Target Detection Complete.');

%% =========================================================
% 3D RECONSTRUCTION
%% =========================================================

disp('Generating 3D Reconstruction...');

point_cloud = reconstruct_3d( ...
                    detections, ...
                    beam_config, ...
                    config);

disp('3D Reconstruction Complete.');

%% =========================================================
% VISUALIZATION ENGINE
%% =========================================================

disp('Launching Visualization Engine...');

render_visualization( ...
        processed_data, ...
        point_cloud, ...
        environment, ...
        config);

disp('Visualization Complete.');

%% =========================================================
% USER INTERFACE
%% =========================================================

disp('Launching Dashboard...');

launch_dashboard( ...
        processed_data, ...
        detections, ...
        point_cloud, ...
        config);

disp('Dashboard Ready.');

%% =========================================================
% HARDWARE INTERFACE
%% =========================================================

disp('Preparing Hardware Interface...');

hardware_interface(config);

disp('Hardware Interface Ready.');

%% =========================================================
% SAVE RESULTS
%% =========================================================

disp('Saving Project Data...');

save_results( ...
        processed_data, ...
        detections, ...
        point_cloud, ...
        config);

disp('Results Saved.');

%% =========================================================
% EXECUTION COMPLETE
%% =========================================================

execution_time = toc(total_runtime);

disp('=================================================');
disp(' SONAR SYSTEM EXECUTION COMPLETE ');
disp('=================================================');

fprintf('Total Runtime : %.2f seconds\n', execution_time);

%% =========================================================
% END OF FILE
%% =========================================================