function [P_db, theta_grid, doa_est, comp_time] = spatial_spectrum(X, array, config, method, varargin)
%SPATIAL_SPECTRUM Unified angle-swept DOA power spectrum, any method.
%
% -------------------------------------------------------------------
% PURPOSE
%   A single, consistent entry point that computes a normalized
%   power-vs-angle spectrum for any of the four beamformers in this
%   module, so callers (in particular COMPARE_BEAMFORMERS.M) do not
%   need to remember four different function signatures just to get
%   "spectrum in dB, DOA estimates, done." Also times the computation
%   for the performance-comparison deliverable.
%
%   This function implements the DS/MVDR/Capon power formulas
%   directly (rather than calling DELAY_SUM/MVDR_BEAMFORMER/
%   CAPON_BEAMFORMER internally) so that, for a given method, the
%   covariance/steering-matrix work is done exactly once regardless
%   of how many angles are scanned. For MUSIC, the eigendecomposition
%   and peak-picking are likewise self-contained here. This does mean
%   the covariance/eigendecomposition logic is intentionally
%   duplicated across this file and the four algorithm-specific
%   files -- see the module-level integration notes for why (keeping
%   each of the nine requested files independently self-contained and
%   callable) and the suggested Phase-2 refactor (extract a shared,
%   private covariance/eigendecomposition utility) if that tradeoff
%   should be revisited.
%
% MATHEMATICAL THEORY
%   See DELAY_SUM.M, MVDR_BEAMFORMER.M / CAPON_BEAMFORMER.M, and
%   MUSIC_DOA.M for the full derivation of each method's power
%   formula. Summary:
%       DS    : P(theta) = a(theta)^H R a(theta) / M^2
%       MVDR  : P(theta) = 1 / ( a(theta)^H R^{-1} a(theta) )
%       Capon : identical formula to MVDR (see MVDR_BEAMFORMER.M for
%               why these are the same algorithm)
%       MUSIC : P(theta) = 1 / ( a(theta)^H E_n E_n^H a(theta) )
%
% ALGORITHM
%   1. Form R_hat = XX^H/N once.
%   2. Dispatch on `method`:
%        'DS'            -> use R_hat directly (no inversion needed)
%        'MVDR' / 'Capon'-> diagonally load R_hat, solve R\A once
%        'MUSIC'         -> eigendecompose R_hat, split subspaces
%   3. Vectorize the angle sweep exactly as in CAPON_BEAMFORMER.M /
%      MUSIC_DOA.M (one [M x M]*[M x K] multiply, not a per-angle
%      loop).
%   4. Normalize to a 0 dB peak; peak-pick with the same
%      dependency-free helper used elsewhere in this module.
%
% INPUT PARAMETERS
%   X      : [M x N] complex array snapshot matrix.
%   array  : struct from ARRAY_GEOMETRY.M.
%   config : struct, needs .fc, .sound_speed, optionally .use_gpu.
%   method : one of 'DS', 'MVDR', 'Capon', 'MUSIC' (case-insensitive).
%
%   Name-Value options:
%       'ThetaGrid'       : [1 x K] angle grid, degrees. Default
%                           -90:0.25:90.
%       'DiagonalLoading' : alpha for MVDR/Capon, default 1e-2.
%       'NumSources'      : D for MUSIC (auto-estimated if omitted;
%                           see MUSIC_DOA.M), or a cap on the number
%                           of reported peaks for DS/MVDR/Capon.
%       'MinSeparationDeg': minimum angular separation between
%                           reported peaks.
%
% OUTPUT PARAMETERS
%   P_db       : [1 x K] spectrum in dB, normalized to a 0 dB peak.
%   theta_grid : the angle grid used (echoed back).
%   doa_est    : row vector of estimated DOA(s), degrees.
%   comp_time  : wall-clock seconds for the covariance/spectrum
%                computation (excludes steering-matrix construction
%                time only if it was already cached -- here it
%                includes it, i.e. this is the full "time to produce
%                a spectrum from raw data," used directly by
%                COMPARE_BEAMFORMERS.M's timing comparison).
%
% COMPLEXITY
%   DS    : O(M^2*N + M^2*K)
%   MVDR/Capon : O(M^2*N + M^3 + M^2*K)
%   MUSIC : O(M^2*N + M^3 + M^2*K)   (eigendecomposition dominates
%           the O(M^3) term in practice)
%
% ASSUMPTIONS
%   Same narrowband / N>>M assumptions as the four algorithm-specific
%   files.
%
% REFERENCES
%   See DELAY_SUM.M, MVDR_BEAMFORMER.M, CAPON_BEAMFORMER.M,
%   MUSIC_DOA.M.
%
% See also DELAY_SUM, MVDR_BEAMFORMER, CAPON_BEAMFORMER, MUSIC_DOA
% -------------------------------------------------------------------

    if nargin < 4
        error('spatial_spectrum:NotEnoughInputs', ...
            'Usage: [P_db,theta_grid,doa_est,comp_time] = spatial_spectrum(X, array, config, method, ...)');
    end

    valid_methods = {'ds','mvdr','capon','music'};
    method_lc = lower(char(method));
    if ~ismember(method_lc, valid_methods)
        error('spatial_spectrum:UnknownMethod', ...
            'Unknown method "%s". Supported: DS, MVDR, Capon, MUSIC.', method);
    end

    p = inputParser;
    addParameter(p, 'ThetaGrid',        -90:0.25:90, @(v) isnumeric(v) && isvector(v));
    addParameter(p, 'DiagonalLoading',  1e-2, @(v) isnumeric(v) && isscalar(v) && v >= 0);
    addParameter(p, 'NumSources',       [],   @(v) isempty(v) || (isnumeric(v) && isscalar(v)));
    addParameter(p, 'MinSeparationDeg', [],   @(v) isempty(v) || (isnumeric(v) && isscalar(v)));
    parse(p, varargin{:});
    theta_grid  = p.Results.ThetaGrid(:)';
    alpha       = p.Results.DiagonalLoading;
    num_sources = p.Results.NumSources;

    validateattributes(X, {'numeric'}, {'2d'}, 'spatial_spectrum', 'X');
    M = size(X,1);
    N = size(X,2);
    if M ~= array.num_elements
        error('spatial_spectrum:DimensionMismatch', ...
            'X has %d rows (channels) but array has %d elements.', M, array.num_elements);
    end

    if isempty(p.Results.MinSeparationDeg)
        min_sep = max(2*mean(diff(theta_grid)), 1);
    else
        min_sep = p.Results.MinSeparationDeg;
    end

    t_start = tic;

    A = steering_vector(array.positions, theta_grid, config.fc, config.sound_speed); % [M x K]
    R_hat = (X * X') / N;

    switch method_lc
        case 'ds'
            U = R_hat * A;                         % [M x K]
            P = real(sum(conj(A) .* U, 1)) / (M^2);

        case {'mvdr','capon'}
            loading = alpha * real(trace(R_hat)) / M;
            R = R_hat + loading * eye(M);
            U = R \ A;
            denom = real(sum(conj(A) .* U, 1));
            P = 1 ./ max(denom, eps);

        case 'music'
            [V, D] = eig(R_hat);
            [eigvals_sorted, idx] = sort(real(diag(D)), 'descend');
            V = V(:, idx);
            if isempty(num_sources)
                D_est = estimate_num_sources(eigvals_sorted, M);
            else
                D_est = num_sources;
                if D_est >= M
                    error('spatial_spectrum:TooManySources', ...
                        'NumSources (%d) must be less than the number of array elements (%d).', D_est, M);
                end
            end
            En = V(:, D_est+1:end);
            Pn = En * En';
            U = Pn * A;
            P = 1 ./ max(real(sum(conj(A) .* U, 1)), eps);
            num_sources = D_est; % for peak-count capping below
    end

    P_db = 10*log10(P / max(P));
    comp_time = toc(t_start);

    doa_est = find_spectral_peaks(P_db, theta_grid, min_sep, num_sources);

end

% -----------------------------------------------------------------------
function D = estimate_num_sources(eigvals_sorted, M)
%ESTIMATE_NUM_SOURCES See MUSIC_DOA.M for full documentation of this
%   identical, intentionally duplicated heuristic.
    max_candidates = max(1, floor(M/2));
    ratios = eigvals_sorted(1:max_candidates) ./ ...
             max(eigvals_sorted(2:max_candidates+1), eps);
    [~, D] = max(ratios);
    D = max(1, min(D, M-1));
end

% -----------------------------------------------------------------------
function locs_deg = find_spectral_peaks(P_db, grid, min_sep_deg, num_sources)
%FIND_SPECTRAL_PEAKS See CAPON_BEAMFORMER.M for full documentation of
%   this identical, intentionally duplicated dependency-free helper.
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
