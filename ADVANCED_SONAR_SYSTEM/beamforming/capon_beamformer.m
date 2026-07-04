function [P_db, theta_grid, doa_est, info] = capon_beamformer(X, array, theta_grid, config, varargin)
%CAPON_BEAMFORMER Capon minimum-variance spatial power spectrum.
%
% -------------------------------------------------------------------
% PURPOSE
%   Computes Capon's (1969) minimum-variance spatial power spectrum
%   over a grid of angles -- a high-resolution alternative to the
%   conventional (Bartlett) spectrum, used here for DOA estimation
%   and visualization rather than for extracting a beamformed signal.
%
% MATHEMATICAL THEORY
%   For each candidate angle theta, Capon's method asks: "what is the
%   minimum output power achievable by a distortionless-at-theta
%   beamformer?" The answer, found by solving the same constrained
%   optimization as MVDR_BEAMFORMER.M (see that file's header for the
%   derivation) and substituting the optimal weight back into
%   w^H R w, is
%
%       P_Capon(theta) = 1 / ( a(theta)^H R^{-1} a(theta) )        (1)
%
%   This is EXACTLY the denominator used to normalize the MVDR weight
%   vector in MVDR_BEAMFORMER.M -- Eq.(1) here and Eq.(2) there are
%   the same mathematical object, just consumed differently (a scalar
%   power spectrum here vs. a signal-extraction weight there). See
%   the note in MVDR_BEAMFORMER.M for the full explanation of why
%   "MVDR" and "Capon" are the same algorithm under two names.
%
%   Resolution: Capon's spectrum is markedly narrower than the
%   conventional (Bartlett) beamformer's mainlobe when sources are
%   well-separated and R is well-estimated, because it adaptively
%   suppresses energy arriving from every angle other than theta
%   when forming the denominator, rather than using fixed weights.
%
% ALGORITHM
%   1. Form sample covariance R_hat = XX^H/N; apply diagonal loading
%      (identical formula to MVDR_BEAMFORMER.M).
%   2. Precompute R^{-1} ONCE (not per angle).
%   3. Build the full steering matrix A = [a(theta_1) ... a(theta_K)]
%      for the whole angle grid in one call to STEERING_VECTOR.M.
%   4. Compute a(theta_k)^H R^{-1} a(theta_k) for ALL K angles at
%      once via U = R^{-1}*A (one [M x M]*[M x K] multiply) followed
%      by a column-wise inner product -- avoiding a per-angle linear
%      solve, which would cost O(M^3*K) instead of O(M^3 + M^2*K).
%   5. P(theta) = 1 ./ diag(A^H * U); convert to dB, normalize to a
%      0 dB peak.
%   6. Locate DOA candidates as local maxima of P_db (see
%      module-level note on the shared, dependency-free peak-picking
%      logic used here and in MUSIC_DOA.M / SPATIAL_SPECTRUM.M).
%
% INPUT PARAMETERS
%   X          : [M x N] complex array snapshot matrix.
%   array      : struct from ARRAY_GEOMETRY.M.
%   theta_grid : [1 x K] angle grid in degrees to scan, e.g.
%                -90:0.25:90. Finer grids cost linearly more time and
%                memory (O(K)) but do not change P(theta) itself.
%   config     : struct, needs .fc, .sound_speed, optionally .use_gpu.
%
%   Name-Value options:
%       'DiagonalLoading' : alpha, default 1e-2 (see MVDR_BEAMFORMER.M).
%       'NumSources'      : if provided, only the strongest
%                           NumSources peaks are returned in doa_est.
%                           If omitted, all peaks clearing a
%                           relative-height floor are returned.
%       'MinSeparationDeg': minimum angular separation enforced
%                           between reported peaks, default 2x the
%                           grid step or 1 degree, whichever is
%                           larger.
%
% OUTPUT PARAMETERS
%   P_db       : [1 x K] Capon spectrum in dB, normalized to a 0 dB
%                peak.
%   theta_grid : the angle grid used (echoed back for convenience).
%   doa_est    : row vector of estimated DOA(s) in degrees.
%   info       : struct with
%                 .R                  loaded covariance used
%                 .diagonal_loading   loading value applied
%
% COMPLEXITY
%   Time  : O(M^2*N) for the covariance, O(M^3) for one matrix solve,
%           O(M^2*K) for the angle sweep. Dominated by O(M^2*K) for
%           fine grids (K >> M), or O(M^2*N) for long snapshot
%           records -- both linear in their respective dimension.
%   Memory: O(M^2) for R plus O(M*K) for the steering matrix.
%
% ASSUMPTIONS
%   Same as MVDR_BEAMFORMER.M (narrowband, N >> M recommended,
%   diagonal loading is a heuristic, not a substitute for adequate
%   snapshot support).
%
% REFERENCES
%   [1] J. Capon, "High-Resolution Frequency-Wavenumber Spectrum
%       Analysis," Proc. IEEE, 57(8):1408-1418, 1969.
%   [2] P. Stoica & R. Moses, "Spectral Analysis of Signals,"
%       Prentice Hall, 2005, Ch. 5.
%
% See also MVDR_BEAMFORMER, MUSIC_DOA, SPATIAL_SPECTRUM
% -------------------------------------------------------------------

    if nargin < 4
        error('capon_beamformer:NotEnoughInputs', ...
            'Usage: [P_db,theta_grid,doa_est,info] = capon_beamformer(X, array, theta_grid, config, ...)');
    end

    p = inputParser;
    addParameter(p, 'DiagonalLoading',  1e-2, @(v) isnumeric(v) && isscalar(v) && v >= 0);
    addParameter(p, 'NumSources',       [],   @(v) isempty(v) || (isnumeric(v) && isscalar(v)));
    addParameter(p, 'MinSeparationDeg', [],   @(v) isempty(v) || (isnumeric(v) && isscalar(v)));
    parse(p, varargin{:});
    alpha       = p.Results.DiagonalLoading;
    num_sources = p.Results.NumSources;

    validateattributes(X, {'numeric'}, {'2d'}, 'capon_beamformer', 'X');
    M = size(X,1);
    N = size(X,2);
    if M ~= array.num_elements
        error('capon_beamformer:DimensionMismatch', ...
            'X has %d rows (channels) but array has %d elements.', M, array.num_elements);
    end

    theta_grid = theta_grid(:)';
    K = numel(theta_grid);

    if isempty(p.Results.MinSeparationDeg)
        grid_step = mean(diff(theta_grid));
        min_sep = max(2*grid_step, 1);
    else
        min_sep = p.Results.MinSeparationDeg;
    end

    use_gpu = isfield(config,'use_gpu') && config.use_gpu && gpu_available();
    if use_gpu
        try
            X = gpuArray(X); %#ok<UNRCH>
        catch
            use_gpu = false;
            warning('capon_beamformer:GPUUnavailable', ...
                'config.use_gpu is true but no usable GPU was found. Falling back to CPU.');
        end
    end

    % ---- covariance + diagonal loading (identical to MVDR_BEAMFORMER.M) ----
    R_hat = (X * X') / N;
    loading = alpha * real(trace(R_hat)) / M;
    R = R_hat + loading * eye(M);

    % ---- steering matrix for the whole grid, one call (vectorized) ----
    A = steering_vector(array.positions, theta_grid, config.fc, config.sound_speed); % [M x K]
    if use_gpu
        A = gpuArray(A); %#ok<UNRCH>
    end

    % ---- Capon spectrum, Eq.(1), vectorized over all K angles ----
    U = R \ A;                                    % [M x K], R^{-1}*A via one solve
    denom = real(sum(conj(A) .* U, 1));            % [1 x K] = a^H R^{-1} a for each angle
    denom = max(denom, eps);                        % guard against numerical negatives
    P = 1 ./ denom;

    if use_gpu
        P = gather(P); %#ok<UNRCH>
    end

    P_db = 10*log10(P / max(P));

    % ---- peak picking (dependency-free; see module notes) ----
    doa_est = find_spectral_peaks(P_db, theta_grid, min_sep, num_sources);

    info = struct();
    info.R = gather_if_needed(R);
    info.diagonal_loading = loading;

end

% -----------------------------------------------------------------------
function locs_deg = find_spectral_peaks(P_db, grid, min_sep_deg, num_sources)
%FIND_SPECTRAL_PEAKS Dependency-free local-maximum peak picker.
%   Deliberately does not use Signal Processing Toolbox's FINDPEAKS
%   so that (a) this module has one less external dependency for
%   its most commonly-called post-processing step, and (b) behavior
%   is identical across MATLAB and Octave (Octave's FINDPEAKS lacks
%   'MinPeakProminence' and rejects negative 'MinPeakHeight', both of
%   which a normalized dB spectrum needs). Duplicated identically in
%   MUSIC_DOA.M and SPATIAL_SPECTRUM.M -- see module-level notes on
%   this deliberate tradeoff.
%
%   Algorithm: find all strict local maxima, sort by height
%   descending, and greedily accept peaks that are at least
%   min_sep_deg away (in angle) from every already-accepted, higher
%   peak. If num_sources is given, return at most that many peaks;
%   otherwise return all peaks within 10 dB of the global maximum
%   (a simple, transparent floor that avoids reporting noise-floor
%   ripple as spurious detections).
    N = numel(P_db);
    if N < 3
        locs_deg = [];
        return;
    end
    is_max = [false, P_db(2:end-1) > P_db(1:end-2) & P_db(2:end-1) >= P_db(3:end), false];
    cand_idx = find(is_max);
    if isempty(cand_idx)
        [~, cand_idx] = max(P_db);
    end

    [~, order] = sort(P_db(cand_idx), 'descend');
    cand_idx = cand_idx(order);

    kept = [];
    for ii = 1:numel(cand_idx)
        idx = cand_idx(ii);
        if isempty(kept) || all(abs(grid(idx) - grid(kept)) >= min_sep_deg)
            if isempty(num_sources) && P_db(idx) < -10
                continue; % below the relative-height floor, skip
            end
            kept(end+1) = idx; %#ok<AGROW>
        end
        if ~isempty(num_sources) && numel(kept) >= num_sources
            break;
        end
    end
    kept = sort(kept);
    locs_deg = grid(kept);
end

% -----------------------------------------------------------------------
function tf = gpu_available()
%GPU_AVAILABLE Local helper: see delay_sum.m for full documentation.
    tf = false;
    if exist('gpuDeviceCount', 'file') == 2 || exist('gpuDeviceCount', 'builtin')
        try
            tf = gpuDeviceCount('available') > 0;
        catch
            tf = false;
        end
    end
end

% -----------------------------------------------------------------------
function out = gather_if_needed(in)
%GATHER_IF_NEEDED Local helper: gather() a gpuArray if applicable,
%   otherwise pass through unchanged. Avoids requiring gpuArray to
%   exist just to check "is this a gpuArray" on a CPU-only machine.
    if isa(in, 'gpuArray')
        out = gather(in); %#ok<UNRCH>
    else
        out = in;
    end
end
