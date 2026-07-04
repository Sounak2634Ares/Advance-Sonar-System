function [P_db, theta_grid, doa_est, eigenvalues, info] = music_doa(X, array, theta_grid, config, varargin)
%MUSIC_DOA MUltiple SIgnal Classification direction-of-arrival estimator.
%
% -------------------------------------------------------------------
% PURPOSE
%   Subspace-based, super-resolution DOA estimator. Exploits the
%   eigenstructure of the array covariance matrix to separate a
%   "signal subspace" (spanned by the source steering vectors) from
%   an orthogonal "noise subspace," giving angular resolution that
%   is not limited by the array's physical aperture the way
%   DELAY_SUM's beamwidth is.
%
% MATHEMATICAL THEORY
%   Data model: X = A(theta) S + N, with A(theta) = [a(theta_1) ...
%   a(theta_D)] the [M x D] matrix of D source steering vectors, S
%   the [D x N] source signals, N spatially-white sensor noise with
%   power sigma^2.
%
%   Covariance eigendecomposition:
%       R = E[xx^H] = A P A^H + sigma^2 I
%   where P = E[ss^H] (source covariance). For D < M sources, R has
%   D "large" eigenvalues (signal-plus-noise) and M-D eigenvalues
%   all equal to sigma^2 (noise-only) in the ideal case. Order the
%   eigenvectors by decreasing eigenvalue:
%       R = V Lambda V^H,  Lambda = diag(lambda_1 >= ... >= lambda_M)
%   The noise subspace is
%       E_n = [ v_{D+1}, ..., v_M ]                                (1)
%   and, crucially, the true source steering vectors are exactly
%   orthogonal to it:  a(theta_i)^H E_n = 0  for i = 1..D.
%
%   MUSIC pseudo-spectrum:
%       P_MUSIC(theta) = 1 / ( a(theta)^H E_n E_n^H a(theta) )     (2)
%   Eq.(2) spikes sharply at the true DOAs because the denominator
%   vanishes there (in the noiseless, infinite-snapshot limit) --
%   this is why MUSIC can resolve sources closer together than the
%   array's conventional (Bartlett) beamwidth: it is not a "beam"
%   being steered at all, but a subspace-orthogonality test.
%
% ALGORITHM
%   1. R_hat = X X^H / N.
%   2. Eigendecompose R_hat; sort eigenvalues descending.
%   3. Split into signal subspace (top D) and noise subspace
%      (remaining M-D) per Eq.(1). D = NumSources if given, else
%      estimated automatically (see below).
%   4. Form the noise-subspace projector Pn = E_n E_n^H once.
%   5. Build the steering matrix for the whole angle grid in one
%      call to STEERING_VECTOR.M; compute a(theta)^H Pn a(theta) for
%      all K angles at once via U = Pn*A followed by a column-wise
%      inner product (identical vectorization strategy to
%      CAPON_BEAMFORMER.M -- one [M x M]*[M x K] multiply instead of
%      a per-angle O(M^2) loop).
%   6. P(theta) = 1 ./ that quantity; convert to dB; peak-pick.
%
%   Automatic source-count estimation (used only if NumSources is
%   not supplied): a simple, transparent eigenvalue-gap heuristic --
%   pick D as the index of the largest ratio drop
%   lambda_i / lambda_{i+1} among the first floor(M/2) eigenvalues.
%   This is a practical default, not a substitute for a proper
%   information-theoretic criterion (AIC/MDL) -- see Known
%   Limitations in the module-level notes.
%
% INPUT PARAMETERS
%   X          : [M x N] complex array snapshot matrix.
%   array      : struct from ARRAY_GEOMETRY.M.
%   theta_grid : [1 x K] angle grid in degrees to scan.
%   config     : struct, needs .fc, .sound_speed, optionally .use_gpu.
%
%   Name-Value options:
%       'NumSources'      : D, the assumed number of sources. If
%                           omitted, estimated automatically (see
%                           above) and reported in `info.num_sources`.
%       'MinSeparationDeg': see CAPON_BEAMFORMER.M; same default rule.
%
% OUTPUT PARAMETERS
%   P_db        : [1 x K] MUSIC pseudo-spectrum in dB, normalized to
%                 a 0 dB peak. NOTE: unlike DS/Capon, MUSIC's P(theta)
%                 is not a physical power in Watts and its absolute
%                 scale carries no meaning -- only peak LOCATIONS are
%                 physically meaningful.
%   theta_grid  : the angle grid used (echoed back).
%   doa_est     : row vector of estimated DOA(s) in degrees, length
%                 D (or fewer, if the grid has fewer resolvable peaks
%                 than D).
%   eigenvalues : [M x 1] eigenvalues of R_hat, sorted descending
%                 (useful for diagnosing the signal/noise subspace
%                 split visually -- a clean split shows a sharp knee).
%   info        : struct with
%                  .num_sources        D actually used
%                  .num_sources_auto   true if D was auto-estimated
%
% COMPLEXITY
%   Time  : O(M^2*N) for the covariance, O(M^3) for the
%           eigendecomposition (dominant cost for large M), O(M^2*K)
%           for the angle sweep.
%   Memory: O(M^2) for R and its eigenvectors, O(M*K) for the
%           steering matrix.
%
% ASSUMPTIONS
%   - Narrowband model (see STEERING_VECTOR.M).
%   - Sources are uncorrelated (or at least not perfectly coherent);
%     fully coherent sources (e.g. specular multipath of the same
%     emitter) collapse the signal subspace rank and MUSIC will
%     under-count sources unless spatial smoothing (not implemented
%     here) is applied first.
%   - D < M (there must be more sensors than sources for a non-empty
%     noise subspace).
%   - N >> M recommended for a reliable eigenvalue split.
%
% REFERENCES
%   [1] R. O. Schmidt, "Multiple Emitter Location and Signal
%       Parameter Estimation," IEEE Trans. Antennas Propag.,
%       34(3):276-280, 1986.
%   [2] H. L. Van Trees, "Optimum Array Processing," Wiley, 2002,
%       Ch. 9.
%
% See also CAPON_BEAMFORMER, MVDR_BEAMFORMER, SPATIAL_SPECTRUM
% -------------------------------------------------------------------

    if nargin < 4
        error('music_doa:NotEnoughInputs', ...
            'Usage: [P_db,theta_grid,doa_est,eigenvalues,info] = music_doa(X, array, theta_grid, config, ...)');
    end

    p = inputParser;
    addParameter(p, 'NumSources',       [], @(v) isempty(v) || (isnumeric(v) && isscalar(v) && v >= 1));
    addParameter(p, 'MinSeparationDeg', [], @(v) isempty(v) || (isnumeric(v) && isscalar(v)));
    parse(p, varargin{:});
    num_sources = p.Results.NumSources;

    validateattributes(X, {'numeric'}, {'2d'}, 'music_doa', 'X');
    M = size(X,1);
    N = size(X,2);
    if M ~= array.num_elements
        error('music_doa:DimensionMismatch', ...
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
            warning('music_doa:GPUUnavailable', ...
                'config.use_gpu is true but no usable GPU was found. Falling back to CPU.');
        end
    end

    % ---- covariance + eigendecomposition ----
    R_hat = (X * X') / N;
    [V, D] = eig(R_hat);
    [eigvals_sorted, idx] = sort(real(diag(D)), 'descend');
    V = V(:, idx);

    num_sources_auto = isempty(num_sources);
    if num_sources_auto
        num_sources = estimate_num_sources(eigvals_sorted, M);
    else
        if num_sources >= M
            error('music_doa:TooManySources', ...
                'NumSources (%d) must be less than the number of array elements (%d).', ...
                num_sources, M);
        end
    end

    En = V(:, num_sources+1:end);    % [M x (M-D)] noise subspace
    if use_gpu
        En = gpuArray(En); %#ok<UNRCH>
    end
    Pn = En * En';                    % [M x M] noise-subspace projector

    % ---- steering matrix for the whole grid, one call (vectorized) ----
    A = steering_vector(array.positions, theta_grid, config.fc, config.sound_speed); % [M x K]
    if use_gpu
        A = gpuArray(A); %#ok<UNRCH>
    end

    % ---- MUSIC pseudo-spectrum, Eq.(2), vectorized over all K angles ----
    U = Pn * A;                                    % [M x K]
    denom = real(sum(conj(A) .* U, 1));             % [1 x K]
    denom = max(denom, eps);
    P = 1 ./ denom;

    if use_gpu
        P = gather(P); %#ok<UNRCH>
        eigvals_sorted = gather(eigvals_sorted); %#ok<UNRCH>
    end

    P_db = 10*log10(P / max(P));

    doa_est = find_spectral_peaks(P_db, theta_grid, min_sep, num_sources);

    eigenvalues = eigvals_sorted;
    info = struct();
    info.num_sources      = num_sources;
    info.num_sources_auto = num_sources_auto;

end

% -----------------------------------------------------------------------
function D = estimate_num_sources(eigvals_sorted, M)
%ESTIMATE_NUM_SOURCES Simple eigenvalue-gap heuristic for source count.
%   Picks D as the index of the largest ratio drop
%   lambda_i / lambda_{i+1} among the first floor(M/2) eigenvalues.
%   A pragmatic default -- not a substitute for AIC/MDL; documented
%   as a Known Limitation and easily overridden via 'NumSources'.
    max_candidates = max(1, floor(M/2));
    ratios = eigvals_sorted(1:max_candidates) ./ ...
             max(eigvals_sorted(2:max_candidates+1), eps);
    [~, D] = max(ratios);
    D = max(1, min(D, M-1));
end

% -----------------------------------------------------------------------
function locs_deg = find_spectral_peaks(P_db, grid, min_sep_deg, num_sources)
%FIND_SPECTRAL_PEAKS Dependency-free local-maximum peak picker.
%   See CAPON_BEAMFORMER.M for full documentation of this
%   intentionally-duplicated helper (identical logic in both files
%   and in SPATIAL_SPECTRUM.M).
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
                continue;
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
