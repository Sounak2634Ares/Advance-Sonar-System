%% =========================================================
% ADVANCED MULTI-BEAM SONAR SYSTEM
% MASTER CONTROL FILE
%
% File        : main.m
% Purpose     : Full Phase 3 Closed-Loop Integration Orchestrator
% Author      : Sounak Jana
% MATLAB      : R2024b
%
% Description :
% Controls the entire sonar execution framework. Tracks operational
% flow, enforces module contracts, and passes real-time matrix telemetry
% from signal conditioning up to the 3D visualization dashboard.
%
% Governing assumptions (held consistent with Phases 1-3a):
%   - Narrowband array approximation: fractional bandwidth BW/fc is
%     small enough that a single steering vector per look angle is valid
%     across the pulse bandwidth.
%   - Far-field / plane-wave incidence at the ULA (no near-field
%     curvature correction).
%   - Baseband analytic (complex envelope) signal representation
%     throughout transmit, propagation, and matched filtering.
%   - Round-trip (monostatic, two-way) delay convention: tau = 2R/c.
%   - Broadside bearing reference (theta = 0 along array normal),
%     symmetric aperture centered at the array midpoint.
%% =========================================================

clc;
clear;
close all;
format compact;

disp('=================================================');
disp(' ADVANCED MULTI-BEAM SONAR SYSTEM — PHASE 3 RUN  ');
disp('=================================================');

%% =========================================================
% PROJECT INTEGRATION PATHS LOADING
%% =========================================================
disp('Loading Project Paths...');

% Only folders that actually exist in this repository are added.
% (Reserved subsystem folders such as 'ui' and 'utils' are added
% conditionally so the framework can absorb them later without
% main.m needing to change.)
required_dirs = { ...
    'config', 'signal_generation', 'environment', 'targets', ...
    'beamforming', 'propagation', 'processing', 'reconstruction', ...
    'visualization', 'gpu', 'hardware', 'data'};

for k = 1:numel(required_dirs)
    if isfolder(required_dirs{k})
        addpath(genpath(required_dirs{k}));
    else
        warning('main:MissingModuleDir', ...
            'Expected module directory "%s" was not found and was skipped.', required_dirs{k});
    end
end

optional_dirs = {'ui', 'utils'};
for k = 1:numel(optional_dirs)
    if isfolder(optional_dirs{k})
        addpath(genpath(optional_dirs{k}));
    end
end

total_runtime = tic;

%% =========================================================
% SYSTEM PARAMETERIZATION & SYNTHESIS SETUP
%% =========================================================
disp('Loading System Configurations...');
% Initialize system variables with uniform linear array parameters
config = struct();
config.Fs = 48e3;                % Sampling rate (Hz)
config.fc = 8e3;                 % Center frequency (Hz)
config.BW = 2e3;                 % Bandwidth (Hz)
config.T  = 0.01;                % Pulse duration (s)
config.N  = 16;                  % Number of array elements
config.c  = 1500;                % Speed of sound (m/s)
% Element spacing fixed at half-wavelength (lambda/2), derived from fc and c
% rather than hardcoded. This is the standard ULA design point: it avoids
% spatial aliasing (grating lobes) while maximizing angular resolution for a
% given element count. The previous d = 0.02 m gave an aperture of only
% (N-1)*d = 0.14 m -- smaller than one wavelength (lambda = c/fc = 0.1875 m)
% at fc = 8 kHz -- producing a ~75 deg mainlobe wide enough to swallow the
% entire CFAR training window and make target detection impossible
% regardless of SNR. With N=16 and d=lambda/2, aperture grows to ~1.4 m and
% the mainlobe narrows to roughly 7 deg, verified numerically to yield real
% CFAR detections at both configured target bearings.
config.d  = config.c / (2 * config.fc); % Element spacing (m), lambda/2
config.SNR_dB = 20;              % Operational SNR target value
config.NumBeams = 64;            % Spatial angular look directions
config.beamformer_type = 'MVDR'; % Choice of spatial processor: 'DAS' or 'MVDR'
config.pfa = 1e-4;               % Constant False Alarm Rate bounds
config.use_gpu = false;          % GPU acceleration toggle (see gpu/initialize_gpu.m)

% Physically-justified TVG absorption coefficient (Thorp's formula,
% evaluated at fc, f in kHz):
%   alpha(f) = 0.11*f^2/(1+f^2) + 44*f^2/(4100+f^2) + 2.75e-4*f^2 + 0.003  [dB/km]
% This replaces signal_processing.m's generic 0.11 dB/km fallback with a
% value derived from the actual operating frequency.
f_kHz = config.fc / 1e3;
config.alpha = 0.11*f_kHz^2/(1+f_kHz^2) + 44*f_kHz^2/(4100+f_kHz^2) ...
             + 2.75e-4*f_kHz^2 + 0.003;   % dB/km

% Spatial (grating-lobe) Nyquist check for the narrowband ULA assumption:
% d <= lambda/2 guarantees an unambiguous visible region over +/-90 deg.
lambda = config.c / config.fc;
if config.d > lambda/2
    warning('main:SpatialAliasingRisk', ...
        ['Element spacing d = %.4f m exceeds lambda/2 = %.4f m at fc = %.0f Hz. ' ...
         'Grating lobes may appear within the steered angular range.'], ...
        config.d, lambda/2, config.fc);
end

disp('Preparing Hardware Interface...');
hardware_interface(config);

disp('Checking GPU Availability...');
try
    gpu_status = initialize_gpu(config); %#ok<NASGU>
catch ME
    warning('main:GPUInitFailed', 'GPU initialization skipped: %s', ME.message);
end

disp('Generating Sonar Waveforms & Target Environment Simulation...');

% -------------------------------------------------------------------
% Transmit reference (matched-filter template) construction.
%
% CRITICAL FIX: this must be the exact complex-conjugate replica of the
% waveform embedded by generate_sonar_echoes.m, i.e. the same analytic
% LFM phase law:
%     phi(t) = 2*pi*( fc*t + (BW/(2T))*t^2 ),   t in [0, T]
%     s(t)   = exp(j*phi(t))
%
% The previous implementation called MATLAB's built-in chirp(), which
% (a) starts its linear-FM sweep at fc - BW/2 instead of fc, introducing
%     a BW/2 frequency offset relative to the actual return model, and
% (b) returns a real-valued cosine rather than the complex analytic
%     signal the receive model and matched filter both expect.
% Both errors reduce matched-filter (pulse-compression) gain and distort
% range sidelobes without raising an error, since dimensions still align.
% -------------------------------------------------------------------
t_chirp = (0:1/config.Fs:config.T).';
tx_signal = exp(1j * 2*pi * (config.fc*t_chirp + (config.BW/(2*config.T))*t_chirp.^2));

% Define target profile arrays
target_list(1) = struct('range', 100, 'bearing', 10,  'target_strength', 1.0); %
target_list(2) = struct('range', 130, 'bearing', -20, 'target_strength', 0.7); %

% Generate raw multi-channel acoustic array returns
[sensor_data, time_vector, params] = generate_sonar_echoes(config, target_list); %#ok<ASGLU>

% Mimic structure packaging from early pipeline setups if required by processing modules
rx_data_struct.signal = sensor_data;

%% =========================================================
% CLOSED-LOOP INTERFACE EXECUTION CHAIN
%% =========================================================
disp('Executing Phase 2 Matrix Signal Conditioners...');
% Runs Fast-Time Matched Filtering and Time-Varying Gain (TVG) normalization
[processed_matrix, proc_meta] = signal_processing(rx_data_struct, tx_signal, config); %#ok<ASGLU>

disp('Executing Phase 3a/b Vectorized Spatial Beamformer...');
% Transforms range-compressed array data from element-space to beam-space
beam_energy = beamforming(processed_matrix, tx_signal, config);

disp('Executing Phase 3b Adaptive 2D CA-CFAR Tracking Search...');
% Extracts valid target returns using adaptive spatial background noise parsing
detections = target_detection(beam_energy, config);

disp('Executing Phase 3c Geographic Target 3D Point-Cloud Mapping...');
% Maps telemetry indices onto explicit 3D space Cartesian geometries
environment_model = struct('Z0', -15); % Establishes a static -15m seabed baseline
pointsXYZ = reconstruct_3d(detections, environment_model);

%% =========================================================
% RESULTS PERSISTENCE
%% =========================================================
disp('Saving Results...');
try
    save_results(processed_matrix, detections, pointsXYZ, config);
catch ME
    warning('main:SaveResultsFailed', 'Could not save results: %s', ME.message);
end

%% =========================================================
% VISUALIZATION
%% =========================================================
% Build axis lookup grids for coordinate conversion sweeps
theta_grid = linspace(-pi/3, pi/3, config.NumBeams);
ranges_grid = (0:size(beam_energy, 1)-1).' * (config.c / (2 * config.Fs));

args_update = struct(...
    'beam_energy', beam_energy, ...
    'theta_vec', theta_grid, ...
    'range_vec', ranges_grid, ...
    'detections', detections); %

% The interactive uifigure-based dashboard requires a graphics-capable
% session (it will error under headless/-batch MATLAB). Guarded so a
% missing display can never prevent the pipeline from completing and
% reporting results; falls back to the static renderer when possible.
try
    disp('Launching Dashboard Graphical Visualization Monitors...');
    app_dashboard = launch_dashboard(environment_model); %
    app_dashboard.runOnceFcn(args_update); %
catch ME
    warning('main:DashboardUnavailable', ...
        'Interactive dashboard unavailable (%s). Falling back to static renderer.', ME.message);
    try
        render_visualization(beam_energy, theta_grid, ranges_grid, pointsXYZ, environment_model, struct());
    catch ME2
        warning('main:StaticRenderFailed', 'Static visualization also failed: %s', ME2.message);
    end
end

%% =========================================================
% REPORT STATISTICS AND PERFORMANCE LOGGING COMPLETION
%% =========================================================
execution_time = toc(total_runtime);
disp('=================================================');
disp(' PHASE 3 COMPLETE — OPERATION NETWORK FULLY ALIGNED ');
disp('=================================================');
fprintf('Active Targets Tracked  : %d positions locked.\n', size(pointsXYZ, 1));
fprintf('Total System Latency    : %.4f seconds\n', execution_time);
