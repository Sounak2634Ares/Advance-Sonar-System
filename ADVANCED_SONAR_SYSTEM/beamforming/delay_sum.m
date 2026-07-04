function [y, w, info] = delay_sum(X, array, theta_look, config)
%DELAY_SUM Conventional (Bartlett) delay-and-sum beamformer.
%
% -------------------------------------------------------------------
% PURPOSE
%   Baseline, non-adaptive beamformer. Steers a fixed, data-
%   independent set of phase weights toward theta_look and sums the
%   array channels. This is the reference every adaptive method
%   (MVDR, Capon, MUSIC) is compared against.
%
% MATHEMATICAL THEORY
%   Weight vector (steered, unity gain at the look direction):
%       w(theta) = a(theta) / M                                   (1)
%   since a(theta)^H a(theta) = M for a unit-modulus steering vector,
%   giving w^H a(theta) = 1 (distortionless at the look direction).
%
%   Beamformer output (time series):
%       y(t) = w^H x(t)  =  (1/M) * sum_{n=1}^{M} conj(a_n) x_n(t)  (2)
%
%   Array (processing) gain for spatially-white noise:
%       G = M   (i.e. 10*log10(M) dB improvement in output SNR over
%       a single element), because coherent signal components sum
%       in amplitude (power grows as M^2) while incoherent noise sums
%       in power (grows as M), giving a net M-fold power SNR gain.
%
%   Delay-and-sum equivalence: see STEERING_VECTOR.M header -- for a
%   narrowband signal, phase-shift weighting (1) is exactly equivalent
%   to applying a time delay tau_n = -(p_n . k_hat)/c per element
%   before summing.
%
% ALGORITHM
%   1. Build the steering matrix for the requested look angle(s).
%   2. Form conventional weights w = A/M (per look angle).
%   3. y = w^H * X   (matrix multiply, vectorized over both channels
%      and, if theta_look has K entries, over all K look directions
%      at once).
%
% INPUT PARAMETERS
%   X          : [M x N] complex array snapshot matrix. Rows = array
%                elements/channels, columns = time snapshots. This
%                row/column convention is used by every function in
%                this module.
%   array      : struct from ARRAY_GEOMETRY.M (needs .positions).
%   theta_look : scalar or [1 x K] vector of look/steering angles in
%                degrees (measured from broadside).
%   config     : struct, needs .fc and .sound_speed.
%
% OUTPUT PARAMETERS
%   y    : [K x N] beamformed output signal(s), one row per requested
%          look angle (K=1 if theta_look is scalar).
%   w    : [M x K] weight vectors used, one column per look angle.
%   info : struct with diagnostic fields:
%           .array_gain_db      theoretical gain 10*log10(M), dB
%           .num_elements       M
%           .num_snapshots      N
%
% COMPLEXITY
%   Time  : O(M*K) to build steering vectors + O(M*N*K) for the
%           matrix multiply w^H*X (dominant term for large N).
%   Memory: O(M*K) for weights + O(K*N) for the output signal(s).
%
% ASSUMPTIONS
%   - Narrowband phase-shift approximation (see STEERING_VECTOR.M).
%   - Spatially white sensor noise (for the stated array-gain figure;
%     the beamformer itself does not require this, only the gain
%     figure quoted in `info` assumes it).
%
% REFERENCES
%   [1] H. L. Van Trees, "Optimum Array Processing," Wiley, 2002,
%       Sec. 2.3 (conventional / Bartlett beamformer).
%
% See also STEERING_VECTOR, MVDR_BEAMFORMER, BEAM_PATTERN
% -------------------------------------------------------------------

    if nargin < 4
        error('delay_sum:NotEnoughInputs', ...
            'Usage: [y,w,info] = delay_sum(X, array, theta_look, config)');
    end

    validateattributes(X, {'numeric'}, {'2d'}, 'delay_sum', 'X');
    M = size(X,1);
    N = size(X,2);

    if M ~= array.num_elements
        error('delay_sum:DimensionMismatch', ...
            'X has %d rows (channels) but array has %d elements.', ...
            M, array.num_elements);
    end

    use_gpu = isfield(config,'use_gpu') && config.use_gpu && gpu_available();
    if use_gpu
        try
            X = gpuArray(X); %#ok<UNRCH> -- exercised only when a GPU exists
        catch
            use_gpu = false;
            warning('delay_sum:GPUUnavailable', ...
                'config.use_gpu is true but no usable GPU was found. Falling back to CPU.');
        end
    end

    A = steering_vector(array.positions, theta_look, config.fc, config.sound_speed); % [M x K]
    if use_gpu
        A = gpuArray(A); %#ok<UNRCH>
    end

    w = A / M;              % [M x K], conventional (Bartlett) weights
    y = w' * X;              % [K x N]

    if use_gpu
        y = gather(y); %#ok<UNRCH>
        w = gather(w); %#ok<UNRCH>
    end

    info = struct();
    info.array_gain_db = 10*log10(M);
    info.num_elements  = M;
    info.num_snapshots = N;

end

% -----------------------------------------------------------------------
function tf = gpu_available()
%GPU_AVAILABLE Local helper: true only if Parallel Computing Toolbox's
%   gpuArray is present AND a usable device can actually be queried.
%   Never throws; always safe to call. Duplicated identically as a
%   local helper in each of the four beamformer files so every file
%   remains independently self-contained (see module-level notes).
    tf = false;
    if exist('gpuDeviceCount', 'file') == 2 || exist('gpuDeviceCount', 'builtin')
        try
            tf = gpuDeviceCount('available') > 0;
        catch
            tf = false;
        end
    end
end
