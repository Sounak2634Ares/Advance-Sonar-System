function [y, w, info] = mvdr_beamformer(X, array, theta_look, config, varargin)
%MVDR_BEAMFORMER Minimum Variance Distortionless Response beamformer.
%
% -------------------------------------------------------------------
% PURPOSE
%   Adaptive (data-dependent) beamformer that minimizes total output
%   power subject to a unity-gain (distortionless) constraint at the
%   look direction. Unlike DELAY_SUM, MVDR uses the measured spatial
%   covariance of the data to place nulls toward interference and
%   noise concentrations while preserving the desired signal.
%
% MATHEMATICAL THEORY
%   Optimization problem:
%       minimize_w   w^H R w
%       subject to   w^H a(theta_look) = 1                        (1)
%   where R = E[x(t) x(t)^H] is the array spatial covariance matrix.
%
%   Solution (Lagrange multipliers / Capon 1969):
%       w_MVDR = R^{-1} a(theta_look) / ( a(theta_look)^H R^{-1} a(theta_look) )   (2)
%
%   NOTE ON NAMING: Eq.(2) is IDENTICAL to the weight vector used by
%   Capon's minimum-variance method (see CAPON_BEAMFORMER.M). "MVDR"
%   and "Capon" are two names, from two literatures, for the same
%   derivation. This file and CAPON_BEAMFORMER.M are kept as separate
%   files (per the requested module structure) but are differentiated
%   by OUTPUT, not by formula:
%       - this file (MVDR_BEAMFORMER) returns the adaptively
%         beamformed TIME-DOMAIN SIGNAL for a chosen look direction
%         (the classical "adaptive beamforming" / signal-extraction
%         use case), while
%       - CAPON_BEAMFORMER.M sweeps Eq.(2) over an angle grid and
%         returns the resulting SPATIAL POWER SPECTRUM (the classical
%         "high-resolution DOA estimation" use case Capon's method is
%         usually taught for).
%   Both are documented here so the equivalence is never hidden.
%
%   Sample covariance & diagonal loading:
%       R_hat = (1/N) * X * X^H                                    (3)
%   R_hat can be poorly conditioned or singular when N < M or when
%   sources are highly correlated/coherent. Diagonal loading
%       R_loaded = R_hat + eps_dl * I,   eps_dl = alpha * trace(R_hat)/M  (4)
%   is applied for numerical stability and robustness to steering-
%   vector/covariance mismatch (a standard, well-established
%   technique -- see Li, Stoica & Wang, 2003).
%
% ALGORITHM
%   1. Form R_hat via (3); apply diagonal loading (4).
%   2. Solve R_loaded * u = a(theta_look) for u  (linear solve, NOT an
%      explicit matrix inverse -- see Computational Complexity).
%   3. w = u / (a(theta_look)^H u).
%   4. y = w^H X.
%
% INPUT PARAMETERS
%   X          : [M x N] complex array snapshot matrix (rows =
%                channels, columns = snapshots).
%   array      : struct from ARRAY_GEOMETRY.M.
%   theta_look : scalar or [1 x K] look angle(s) in degrees.
%   config     : struct, needs .fc, .sound_speed, optionally .use_gpu.
%
%   Name-Value options:
%       'DiagonalLoading' : alpha in Eq.(4), default 1e-2 (1% of the
%                           average eigenvalue of R_hat). Increase
%                           for more robustness / less aggressive
%                           nulling when snapshot count is low or
%                           steering-vector mismatch is a concern;
%                           decrease for maximum interference
%                           rejection when R_hat is well estimated.
%
% OUTPUT PARAMETERS
%   y    : [K x N] adaptively beamformed output signal(s).
%   w    : [M x K] MVDR weight vectors used.
%   info : struct with
%           .R                  [M x M] loaded covariance matrix used
%           .diagonal_loading   the numeric loading value applied
%           .condition_number   cond(R_loaded) (numerical health check)
%
% COMPLEXITY
%   Time  : O(M^2*N) to form R_hat, O(M^3) to solve the linear system
%           per look angle (Cholesky-based backslash), so O(M^3*K)
%           total for K look angles. For typical sonar arrays
%           (M ~ 10-100) this is negligible next to the O(M^2*N) data
%           term when N is large (N here can be ~10^3-10^4 snapshots
%           per pulse at fs = 1 MHz).
%   Memory: O(M^2) for R plus O(M*N) for X.
%
% ASSUMPTIONS
%   - Narrowband phase-shift model (see STEERING_VECTOR.M).
%   - Snapshots in X are drawn from a (locally) wide-sense-stationary
%     process so that the sample covariance (3) is a meaningful
%     estimate of R; N >> M is recommended (rule of thumb N >= 2M,
%     ideally N >= 5-10 M) for a well-conditioned estimate.
%   - Diagonal loading (4) is a heuristic robustness measure, not a
%     substitute for adequate snapshot support.
%
% REFERENCES
%   [1] J. Capon, "High-Resolution Frequency-Wavenumber Spectrum
%       Analysis," Proc. IEEE, 57(8), 1969.
%   [2] H. L. Van Trees, "Optimum Array Processing," Wiley, 2002,
%       Sec. 6.2 (MVDR / Capon beamformer).
%   [3] J. Li, P. Stoica, Z. Wang, "On Robust Capon Beamforming and
%       Diagonal Loading," IEEE Trans. Signal Process., 51(7), 2003.
%
% See also CAPON_BEAMFORMER, STEERING_VECTOR, DELAY_SUM
% -------------------------------------------------------------------

    if nargin < 4
        error('mvdr_beamformer:NotEnoughInputs', ...
            'Usage: [y,w,info] = mvdr_beamformer(X, array, theta_look, config, ...)');
    end

    p = inputParser;
    addParameter(p, 'DiagonalLoading', 1e-2, @(v) isnumeric(v) && isscalar(v) && v >= 0);
    parse(p, varargin{:});
    alpha = p.Results.DiagonalLoading;

    validateattributes(X, {'numeric'}, {'2d'}, 'mvdr_beamformer', 'X');
    M = size(X,1);
    N = size(X,2);

    if M ~= array.num_elements
        error('mvdr_beamformer:DimensionMismatch', ...
            'X has %d rows (channels) but array has %d elements.', ...
            M, array.num_elements);
    end
    if N < M
        warning('mvdr_beamformer:FewSnapshots', ...
            ['Only %d snapshots for a %d-element array (N < M). The sample ' ...
             'covariance will be rank-deficient; relying entirely on ' ...
             'diagonal loading for stability. Consider N >= 2*M.'], N, M);
    end

    use_gpu = isfield(config,'use_gpu') && config.use_gpu && gpu_available();
    if use_gpu
        try
            X = gpuArray(X); %#ok<UNRCH>
        catch
            use_gpu = false;
            warning('mvdr_beamformer:GPUUnavailable', ...
                'config.use_gpu is true but no usable GPU was found. Falling back to CPU.');
        end
    end

    % ---- sample covariance + diagonal loading (Eq. 3-4) ----
    R_hat = (X * X') / N;
    loading = alpha * real(trace(R_hat)) / M;
    R = R_hat + loading * eye(M);

    % ---- steering vector(s) at the look direction(s) ----
    A = steering_vector(array.positions, theta_look, config.fc, config.sound_speed); % [M x K]
    if use_gpu
        A = gpuArray(A); %#ok<UNRCH>
    end
    K = size(A,2);

    % ---- solve R*u = a for each look angle (linear solve, not inv()) ----
    U = R \ A;                                  % [M x K]
    denom = real(sum(conj(A) .* U, 1));          % [1 x K], = a^H R^{-1} a
    w = U ./ denom;                               % [M x K], broadcast divide

    y = w' * X;                                    % [K x N]

    if use_gpu
        y = gather(y); w = gather(w); R = gather(R); %#ok<UNRCH>
    end

    info = struct();
    info.R                = R;
    info.diagonal_loading = loading;
    info.condition_number = cond(R);

end

% -----------------------------------------------------------------------
function tf = gpu_available()
%GPU_AVAILABLE Local helper: see delay_sum.m for full documentation of
%   this identical, intentionally duplicated helper.
    tf = false;
    if exist('gpuDeviceCount', 'file') == 2 || exist('gpuDeviceCount', 'builtin')
        try
            tf = gpuDeviceCount('available') > 0;
        catch
            tf = false;
        end
    end
end
