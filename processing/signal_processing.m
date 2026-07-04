function [processed_data, processing_metadata] = signal_processing(rx_data, tx_signal, config)
% SIGNAL_PROCESSING Perform Multi-Channel Matched Filtering and TVG Compensation.
%
% File        : signal_processing.m
% Purpose     : Phase 2 Signal Conditioning Core (Final Pipeline Integrated)
% Author      : Central Coordinator AI
%
% Inputs:
%   rx_data    - Struct containing simulation fields OR raw [NumSamples x N] matrix
%   tx_signal  - Vector representing reference transmitted waveform
%   config     - Struct containing system configuration parameters
%
% Outputs:
%   processed_data      - [NumSamples x N] range-compressed, gain-equalized matrix
%   processing_metadata - Struct returning applied gains and calculation logs

% =========================================================================
% 0. DEFENSIVE STRUCT DETECTION & FIELD UNWRAPPING
% =========================================================================
if isstruct(rx_data)
    % Map to the verified field from the simulation's structural bundle
    if isfield(rx_data, 'signal')
        rx_matrix = rx_data.signal;
    elseif isfield(rx_data, 'sensor_data')
        rx_matrix = rx_data.sensor_data;
    elseif isfield(rx_data, 'data')
        rx_matrix = rx_data.data;
    else
        fields = fieldnames(rx_data);
        error('SignalProcessing:InvalidStruct', ...
            'rx_data was passed as a struct, but no valid data field was found. Available fields: %s', ...
            strjoin(fields, ', '));
    end
else
    % rx_data is already a standard numeric matrix
    rx_matrix = rx_data;
end

% Ensure arrays are evaluated correctly using the extracted matrix
[NumSamples, N] = size(rx_matrix);

% =========================================================================
% 1. CASE-INSENSITIVE PARAMETER EXTRACTION GATEWAY
% =========================================================================
if isfield(config, 'Fs'),          Fs = config.Fs;
elseif isfield(config, 'fs'),      Fs = config.fs;
elseif isfield(config, 'SampleRate'), Fs = config.SampleRate;
else, error('SignalProcessing:MissingFs', 'Configuration missing sampling rate (Fs).');
end

if isfield(config, 'c'), c = config.c; else, c = 1500; end
if isfield(config, 'alpha'), alpha = config.alpha; else, alpha = 0.11; end

if isfield(config, 'T'),          pulse_duration = config.T;
elseif isfield(config, 't_pulse'), pulse_duration = config.t_pulse;
else, pulse_duration = 0.01; 
end

time_vector = (0:NumSamples-1).' / Fs; 

% =========================================================================
% 2. TIME-VARYING GAIN (TVG) EQUALIZATION
% =========================================================================
ranges = (c * time_vector) / 2;
ranges(ranges < eps) = eps; 

G_dB = 40 * log10(ranges) + 2 * alpha * (ranges / 1000);

clamp_idx = find(time_vector >= pulse_duration, 1, 'first');
if ~isempty(clamp_idx)
    target_gain_scalar = G_dB(clamp_idx);
    blank_mask = time_vector < pulse_duration;
    G_dB(blank_mask) = target_gain_scalar;
else
    G_dB(:) = 0;
end

gain_vector = 10 .^ (G_dB / 20); 

% Multiply using the raw numeric extracted matrix array
tvg_scaled_data = rx_matrix .* gain_vector;

% =========================================================================
% 3. FFT MATCHED FILTERING / PULSE COMPRESSION
% =========================================================================
N_fft = NumSamples; 

tx_signal = tx_signal(:);
tx_len = length(tx_signal);

replica = zeros(N_fft, 1);
if tx_len <= N_fft
    replica(1:tx_len) = tx_signal;
else
    replica = tx_signal(1:N_fft);
end

if isfield(config, 'window_type') && strcmpi(config.window_type, 'Hamming')
    win = hamming(min(tx_len, N_fft));
    replica(1:length(win)) = replica(1:length(win)) .* win;
end

X_rx = fft(tvg_scaled_data, N_fft, 1);
H_ref = fft(replica, N_fft);

Compressed_Spectrum = X_rx .* conj(H_ref);
processed_data = ifft(Compressed_Spectrum, N_fft, 1);

% =========================================================================
% 4. METADATA DIAGNOSTICS LOGGING
% =========================================================================
processing_metadata.gain_curve_dB = G_dB;
processing_metadata.applied_linear_gains = gain_vector;
processing_metadata.fft_points = N_fft;

end