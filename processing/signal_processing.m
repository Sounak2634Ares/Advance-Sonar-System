function [processed_data, processing_metadata] = signal_processing(rx_data, tx_signal, config)
% SIGNAL_PROCESSING Perform Multi-Channel Matched Filtering and TVG Compensation.
%
% File        : signal_processing.m
% Purpose     : Phase 2 Signal Conditioning Core
% Author      : Central Coordinator AI (Gemini)
%
% Inputs:
%   rx_data    - [NumSamples x N] complex-valued baseband raw acoustic matrix
%   tx_signal  - [NumChirpSamples x 1] vector representing reference transmitted waveform
%   config     - Struct containing fields: Fs, c, alpha, window_type (optional)
%
% Outputs:
%   processed_data      - [NumSamples x N] range-compressed, gain-equalized matrix
%   processing_metadata - Struct returning applied gains and calculation logs

% =========================================================================
% 1. VALIDATION AND PARAMETER EXTRACTION
% =========================================================================
[NumSamples, N] = size(rx_data);
Fs = config.Fs;

if isfield(config, 'c')
    c = config.c;
else
    c = 1500; % Fallback speed of sound in m/s
end

if isfield(config, 'alpha')
    alpha = config.alpha;
else
    alpha = 0.11; % Fallback absorption coefficient in dB/km (Thorpe default)
end

% Create unambiguous time vector
time_vector = (0:NumSamples-1).' / Fs; 

% =========================================================================
% 2. PIECEWISE TIME-VARYING GAIN (TVG) EQUALIZATION
% =========================================================================
% Convert two-way time of flight to one-way range (meters)
ranges = (c * time_vector) / 2;

% Apply a tiny epsilon offset to prevent log10(0) singularities at range = 0
ranges(ranges < eps) = eps; 

% Calculate exact Two-Way Transmission Loss: TL = 2 * [20*log10(R) + alpha * R_km]
% Distributed linearly: G_dB = 40*log10(R) + 2 * alpha * (R / 1000)
G_dB = 40 * log10(ranges) + 2 * alpha * (ranges / 1000);

% Define Piecewise TVG limits (Blanking window & saturation ceiling)
t0 = isfield(config, 'T') * config.T; % Use pulse width as standard blanking window
if isempty(t0) || t0 == 0, t0 = 0.002; end % Fallback to 2ms

% Apply Clamping Limits
G_dB(time_vector < t0) = G_dB(find(time_vector >= t0, 1, 'first')); 

% Convert to linear voltage/pressure multiplier
gain_vector = 10 .^ (G_dB / 20); 

% Apply TVG across all N columns using implicit matrix expansion
tvg_scaled_data = rx_data .* gain_vector;

% =========================================================================
% 3. FFT-BASED FREQUENCY DOMAIN MATCHED FILTERING
% =========================================================================
% Match FFT execution block to the exact dimensions of incoming matrix row count
N_fft = NumSamples; 

% Prepare reference chirp vector
tx_len = length(tx_signal);
replica = zeros(N_fft, 1);
replica(1:tx_len) = tx_signal;

% Optional Frequency-Domain Windowing to control Peak Sidelobe Levels (PSLL)
if isfield(config, 'window_type') && strcmpi(config, 'window_type', 'Hamming')
    win = hamming(tx_len);
    replica(1:tx_len) = replica(1:tx_len) .* win;
end

% Compute fast-time 1D Fourier Transform of input signal channels and template
X_rx = fft(tvg_scaled_data, N_fft, 1);
H_ref = fft(replica, N_fft);

% Apply Conjugate Spectral Multiplication (Matched Filter Identity: H(f) = S*(f))
% Utilizing implicit expansion across all N channels
Compressed_Spectrum = X_rx .* conj(H_ref);

% Revert to Time-Domain via Inverse Fast Fourier Transform
processed_data = ifft(Compressed_Spectrum, N_fft, 1);

% Pack metadata metrics for diagnostics
processing_metadata.gain_curve_dB = G_dB;
processing_metadata.applied_linear_gains = gain_vector;
processing_metadata.fft_points = N_fft;

end