function [sensor_data, time_vector, params] = generate_sonar_echoes(sys_config, target_list)
% GENERATE_SONAR_ECHOES Simulates active acoustic returns for a Uniform Linear Array.
%
% File        : generate_sonar_echoes.m
% Purpose     : Phase 3 Verification - Synthetic Environment Matrix Simulator
% Author      : Sounak Jana / Gemini Head Coordinator
%
% Inputs:
%   sys_config  - Struct containing Fs, fc, BW, T, N, d, c, and SNR_dB parameters
%   target_list - Array of structs containing target profiles: 'range', 'bearing', 'target_strength'
%
% Outputs:
%   sensor_data - [NumSamples x N] complex baseband analytic echo returns matrix
%   time_vector - [NumSamples x 1] time vector coordinates mapping array
%   params      - Compilation metadata struct log wrapper

    % 1. Parse configuration parameters with absolute defensive stability
    if isfield(sys_config, 'Fs'),  Fs = sys_config.Fs;  else, Fs = 48e3; end
    if isfield(sys_config, 'fc'),  fc = sys_config.fc;  else, fc = 8e3;  end
    if isfield(sys_config, 'BW'),  BW = sys_config.BW;  else, BW = 2e3;  end
    if isfield(sys_config, 'T'),   T  = sys_config.T;   else, T  = 0.01; end
    if isfield(sys_config, 'N'),   N  = sys_config.N;   else, N  = 8;    end
    if isfield(sys_config, 'd'),   d  = sys_config.d;   else, d  = 0.02; end
    if isfield(sys_config, 'c'),   c  = sys_config.c;   else, c  = 1500; end
    if isfield(sys_config, 'SNR_dB'), SNR_dB = sys_config.SNR_dB; else, SNR_dB = 20; end

    % 2. Establish time boundaries for total matrix grid space
    % Maximum two-way propagation coverage window benchmarked for deep tracking space
    t_max = 0.3; 
    time_vector = (0:1/Fs:t_max).';
    NumSamples = length(time_vector);
    
    sensor_data = zeros(NumSamples, N);
    M = length(target_list);
    
    if M == 0
        % Return pure uncorrelated ambient white noise floor if scene space is empty
        sensor_data = (randn(NumSamples, N) + 1j*randn(NumSamples, N)) * (10^(-SNR_dB/20) / sqrt(2));
        params = sys_config;
        return;
    end

    % 3. Calculate Spatial-Temporal Geometric Delays for Element Channels
    n_indices = (0:N-1) - (N-1)/2; % Symmetric layout offsets relative to array axis center

    for m = 1:M
        r_target = target_list(m).range;
        
        % Convert bearing parameter to radians gracefully whether input is in degrees or radians
        b_target = target_list(m).bearing;
        if abs(b_target) > 2*pi
            b_target = b_target * pi / 180; % Convert from degrees if explicit values detected
        end
        
        ts_amp = 1.0;
        if isfield(target_list(m), 'target_strength'), ts_amp = target_list(m).target_strength; end

        % Two-Way Monostatic Range Time Propagation Delay (Round Trip Time)
        tau_range = 2 * r_target / c;

        for n = 1:N
            % Inter-element array aperture delay offset based on plane wave arrival angles
            tau_elem = (n_indices(n) * d * sin(b_target)) / c;
            
            % --- FIXED THE TYPO: Correct combined variable summation calculation ---
            tau_total = tau_range + tau_elem; 

            % Evaluate time vectors passing threshold bounds for analytic envelope injection
            t_shifted = time_vector - tau_total;
            mask = (t_shifted >= 0) & (t_shifted <= T);

            if any(mask)
                % Synthesize spatial analytic signal return tracking vectors
                % Computes instantaneous LFM chirp frequency profiles coupled to array element offsets
                t_envelope = t_shifted(mask);
                envelope = ts_amp * rectwin(length(t_envelope)); % Clear rectangular temporal shaping
                
                % Compute chirp transmission frequency arguments matching your master tx parameters
                phase_arg = 2 * pi * (fc * t_envelope + (BW / (2 * T)) * t_envelope.^2);
                
                % Inject complex analytic components into matrix tracking channels
                sensor_data(mask, n) = sensor_data(mask, n) + envelope .* exp(1j * phase_arg);
            end
        end
    end

    % 4. Inject Additive White Gaussian Noise Floor matching configured SNR
    Ps = mean(sum(abs(sensor_data).^2, 1) / NumSamples);
    if Ps > 0
        SNR_lin = 10^(SNR_dB / 10);
        Pn = Ps / SNR_lin;
        sigma_noise = sqrt(Pn / 2);
        noise_matrix = sigma_noise * (randn(NumSamples, N) + 1j * randn(NumSamples, N));
        sensor_data = sensor_data + noise_matrix;
    end

    % Export execution metrics parameters
    params = sys_config;
    params.num_targets = M;
    params.sample_grid_dims = size(sensor_data);
end