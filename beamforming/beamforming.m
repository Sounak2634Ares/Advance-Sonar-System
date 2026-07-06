function beam_energy = beamforming(processed_data, tx_signal, config)
% BEAMFORMING Spatial processing engine master dispatcher module.
%
% File        : beamforming.m
% Purpose     : Phase 3a/b Spatial Array Processing Core Interface
% Author      : Gemini Head Coordinator / ChatGPT Architecture
%
% Inputs:
%   processed_data - [NumSamples x N] numeric complex range-compressed matrix
%   tx_signal      - Vector representing transmitted reference signal template
%   config         - System configuration structure
%
% Outputs:
%   beam_energy    - [NumSamples x NumBeams] numeric energy power matrix

% =========================================================================
% 1. DEFENSIVE RIGOROUS SPECIFICATION GUARD RAILS
% =========================================================================
validateattributes(processed_data, {'single', 'double'}, {'2d', 'nonempty'});

% Dynamic configuration channel mapping (handles 'N' or 'n')
if isfield(config, 'N'),         N_elements = config.N;
elseif isfield(config, 'n'),     N_elements = config.n;
else,                            N_elements = size(processed_data, 2);
end

assert(size(processed_data, 2) == N_elements, ...
    'Beamforming:ElementMismatch', 'Input matrix column count must match configuration element count N.');

assert(all(isfinite(processed_data(:))), ...
    'Beamforming:NumericalInstability', 'Input range-compressed matrix contains non-finite values (NaN/Inf).');

% Establish structural fallback configurations if fields are omitted
if ~isfield(config, 'beamformer_type'), config.beamformer_type = 'DAS'; end
if ~isfield(config, 'NumBeams'),        config.NumBeams = 64; end

% =========================================================================
% 2. EXECUTION STRATEGY DISPATCHER ROUTINE
% =========================================================================
switch upper(string(config.beamformer_type))
    case 'DAS'
        beam_energy = beamforming_DAS(processed_data, config);
        
    case 'MVDR'
        beam_energy = beamforming_MVDR(processed_data, config);
        
    otherwise
        error('Beamforming:UnknownType', ...
            'Execution aborted. Unrecognized beamformer strategy selection: "%s". Use "DAS" or "MVDR".', ...
            config.beamformer_type);
end

end

% =========================================================================
% 3. HIGH-PERFORMANCE VECTORIZED DELAY-AND-SUM ENGINE (DAS)
% =========================================================================
function beam_energy = beamforming_DAS(X, config)
[NumSamples, N] = size(X);
NumBeams = config.NumBeams;

% Extract environmental parameters with defensive defaults
if isfield(config, 'fc'), fc = config.fc; else, fc = 8e3; end
if isfield(config, 'c'),  c = config.c;   else, c = 1500; end
if isfield(config, 'd'),  d = config.d;   else, d = c / (2 * fc); end

% Generate steering angular space over [-60, +60] degrees Field-of-View
theta_vec = linspace(-pi/3, pi/3, NumBeams); 

% Compute element layout offsets relative to geometric center phase-reference
n_vec = (0:N-1).' - (N-1)/2; % [N x 1] column vector

% Vectorized Array Manifold Generation: [N x NumBeams]
% [a(theta)]_n = exp(-j * 2 * pi * fc * (n_offset * d * sin(theta)) / c)
% This matches the physical receive-side phase model in generate_sonar_echoes.m.
steering_matrix = exp(-1j * 2 * pi * fc * (n_vec * d * sin(theta_vec)) / c);

% Compute spatial Bartlett weights: w = (1/N) * a*(theta), i.e. the CONJUGATE
% of the array manifold, so that combining coherently cancels the receive-side
% phase at theta = true bearing. (Using the manifold un-conjugated here mirrors
% the recovered bearing about broadside -- verified empirically.) -> [N x NumBeams]
W_DAS = (1 / N) * conj(steering_matrix);

% Execute High-Speed Fast-Time Matrix Multiplications via optimized BLAS pipeline
% [NumSamples x N] * [N x NumBeams] = [NumSamples x NumBeams]
beam_output = X * W_DAS;

% Compute instantaneous power response magnitude
beam_energy = abs(beam_output).^2;
end

% =========================================================================
% 4. ADAPTIVE REGULARIZED MINIMUM VARIANCE DISTORTIONLESS RESPONSE (MVDR)
% =========================================================================
function beam_energy = beamforming_MVDR(X, config)
[NumSamples, N] = size(X);
NumBeams = config.NumBeams;

if isfield(config, 'fc'), fc = config.fc; else, fc = 8e3; end
if isfield(config, 'c'),  c = config.c;   else, c = 1500; end
if isfield(config, 'd'),  d = config.d;   else, d = c / (2 * fc); end
if isfield(config, 'mvdr_epsilon'), epsilon = config.mvdr_epsilon; else, epsilon = 1e-2; end

theta_vec = linspace(-pi/3, pi/3, NumBeams);
n_vec = (0:N-1).' - (N-1)/2;

% Precompute Steering Vectors: [N x NumBeams]
% Sign convention: beam_energy is formed below as X*w_mvdr (not X*conj(w_mvdr)),
% so for the distortionless constraint a_theta^H*w = 1 to correspond to a peak
% at theta = true bearing, a_theta must be the CONJUGATE of the physical
% receive-side manifold used in generate_sonar_echoes.m (exp(-j*2*pi*fc*n*d*sin(bearing)/c)).
% Verified empirically against a known single-target bearing.
A_matrix = exp(1j * 2 * pi * fc * (n_vec * d * sin(theta_vec)) / c);

% Estimate Global Spatial Sample Covariance Matrix over all fast-time snapshot windows
% R_xx shape: [N x N]
R_xx = (X' * X) / NumSamples;

% Inject Claude's Dynamic Scale-Invariant Diagonal Loading Regularization
trace_val = trace(R_xx);
gamma = epsilon * (trace_val / N);
R_xx_reg = R_xx + gamma * eye(N, 'like', R_xx);

% Preallocate output beam space allocation matrix
beam_energy = zeros(NumSamples, NumBeams, 'like', real(X));

% Process beams. Avoid explicit matrix inversion using highly accurate backslash solver
for b = 1:NumBeams
    a_theta = A_matrix(:, b); % Get current steering vector [N x 1]
    
    % Compute R_inv * a_theta via linear system solving equations
    R_inv_a = R_xx_reg \ a_theta; 
    
    % Compute Capon normalization denominator: a^H * R_inv * a
    denom = a_theta' * R_inv_a;
    
    % Compute normalized spatial adaptive weight optimization vectors
    w_mvdr = R_inv_a / denom;
    
    % Linearly project computed weights across multi-channel matrix streams
    beam_energy(:, b) = abs(X * w_mvdr).^2;
end

end