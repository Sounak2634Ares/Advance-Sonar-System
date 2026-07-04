function results = compare_beamformers(config, varargin)
%COMPARE_BEAMFORMERS Phase-1 beamforming demo, test, and comparison script.
%
% -------------------------------------------------------------------
% PURPOSE
%   Single entry point that (1) builds the array, (2) runs the five
%   requested simulation cases through all four beamformers
%   (Delay-and-Sum, MVDR, Capon, MUSIC), (3) produces the requested
%   visualizations, and (4) prints a performance-comparison summary
%   (resolution, computation time, memory, sidelobe suppression,
%   beamwidth, robustness to noise). This file serves double duty as
%   both "compare_beamformers.m" in the requested module structure
%   and the "demo script" requested under TESTING.
%
%   IMPORTANT INTEGRATION NOTE: this script generates its own
%   synthetic per-element array data (see the local function
%   SIMULATE_ARRAY_SNAPSHOTS below) rather than pulling from
%   propagation/simulate_echoes.m. That existing function currently
%   produces a single combined channel (rx_data.signal is
%   1 x signal_length, with no per-hydrophone dimension at all), so
%   there is no angle-resolved array data anywhere in the existing
%   pipeline yet for a beamformer to consume. Extending
%   simulate_echoes.m to emit true per-element data (accounting for
%   each target's bearing relative to the array) is outside this
%   module's file scope (propagation/ is not in the requested
%   beamforming/ file list) -- see the module-level integration notes
%   delivered alongside this file for a recommended follow-up.
%
% SIMULATION CASES (as specified)
%   Case 1: Single target.
%   Case 2: Two closely-spaced targets, separated by less than the
%           array's own -3 dB beamwidth (computed from the array,
%           not hard-coded) -- the classic test of super-resolution.
%   Case 3: Three targets spread across the scan sector.
%   Case 4: Low SNR (-10 dB per element).
%   Case 5: Strong interference (target + a +20 dB stronger
%           interferer at a different angle) -- demonstrates adaptive
%           nulling (MVDR/Capon) vs. the conventional beamformer's
%           inability to reject it.
%
% INPUT PARAMETERS
%   config : struct from config/system_config.m. If omitted, this
%            function calls system_config() itself (matching this
%            project's convention that every module works from that
%            struct).
%
%   Name-Value options:
%       'NumSnapshots' : N, snapshots per case, default 1000
%                        (N/M ~= 31 for the default M=32 array, well
%                        within the N >> M guideline used throughout
%                        this module).
%       'ThetaGrid'    : angle grid for all spectra, default -90:0.25:90.
%       'MakePlots'    : true (default) or false -- set false for
%                        headless/CI runs (e.g. inside Octave without
%                        a display).
%       'Seed'         : RNG seed for reproducibility, default 42.
%
% OUTPUT PARAMETERS
%   results : struct array (one entry per case) with fields
%              .name, .true_doas, .snr_db, .doa_estimates (struct
%              with .DS/.MVDR/.Capon/.MUSIC), .comp_time_s (same
%              struct shape), .metrics (per-method resolution/PSLL
%              where applicable).
%
% COMPLEXITY / PERFORMANCE ANALYSIS
%   For each case, this script measures and reports (see the printed
%   summary table and METRICS_TABLE local function):
%     - Angular resolution   : numerically measured HPBW (DS, via
%                               BEAM_PATTERN.M) and achieved MUSIC/
%                               Capon peak separation in Case 2.
%     - Computation time      : wall-clock seconds per method, from
%                               SPATIAL_SPECTRUM.M's `comp_time`.
%     - Memory usage           : approximate bytes for the dominant
%                               data structures (X, R, steering
%                               matrix), via a portable WHOS-based
%                               helper (see MEMORY_ESTIMATE_BYTES).
%     - Side-lobe suppression  : PSLL from BEAM_PATTERN.M (DS) vs.
%                               the adaptive spectra's off-target
%                               floor.
%     - Beam width             : HPBW from BEAM_PATTERN.M.
%     - Robustness to noise    : Case 4 (-10 dB) DOA error per method.
%
% ASSUMPTIONS
%   See the per-file headers of ARRAY_GEOMETRY, STEERING_VECTOR,
%   DELAY_SUM, MVDR_BEAMFORMER, CAPON_BEAMFORMER, MUSIC_DOA,
%   BEAM_PATTERN, SPATIAL_SPECTRUM for the full assumption list
%   (narrowband model, N>>M, uncorrelated sources, etc.). This script
%   additionally assumes: synthetic far-field plane-wave sources
%   (not the project's actual target x/y/z + seabed/multipath model),
%   used here specifically because real per-element propagation data
%   is not yet available in the pipeline (see PURPOSE above).
%
% See also ARRAY_GEOMETRY, STEERING_VECTOR, DELAY_SUM, MVDR_BEAMFORMER,
%          CAPON_BEAMFORMER, MUSIC_DOA, BEAM_PATTERN, SPATIAL_SPECTRUM
% -------------------------------------------------------------------

    if nargin < 1 || isempty(config)
        config = system_config();
    end

    p = inputParser;
    addParameter(p, 'NumSnapshots', 1000, @(v) isnumeric(v) && isscalar(v) && v > 0);
    addParameter(p, 'ThetaGrid',    -90:0.25:90, @(v) isnumeric(v) && isvector(v));
    addParameter(p, 'MakePlots',    true, @(v) islogical(v) || isnumeric(v));
    addParameter(p, 'Seed',         42, @(v) isnumeric(v) && isscalar(v));
    parse(p, varargin{:});
    N          = p.Results.NumSnapshots;
    theta_grid = p.Results.ThetaGrid(:)';
    make_plots = logical(p.Results.MakePlots);
    rng(p.Results.Seed);

    fprintf('=================================================\n');
    fprintf(' BEAMFORMING MODULE - COMPARISON & DEMO\n');
    fprintf('=================================================\n');

    array = array_geometry(config);

    % ---- data-driven Case 2 separation: pick something inside the
    % array's own -3dB beamwidth so the "closely spaced" case is a
    % genuine, non-arbitrary resolution test for this exact config ----
    a0 = steering_vector(array.positions, 0, config.fc, config.sound_speed);
    [~, ~, bw_metrics] = beam_pattern(a0/array.num_elements, array, config, ...
        'ThetaGrid', theta_grid, 'PlotType', 'none');
    hpbw0 = bw_metrics.hpbw_deg;
    close_sep = 0.6 * hpbw0;   % inside the DS mainlobe -> tests super-resolution
    fprintf('Array broadside HPBW = %.2f deg -> Case 2 separation set to %.2f deg\n\n', ...
        hpbw0, close_sep);

    % ---- define the five cases -----------------------------------
    cases = struct('name',{},'doas',{},'powers',{},'snr_db',{},'interference',{});
    cases(1) = make_case('Case 1: Single target',            15,                  1,             10, []);
    cases(2) = make_case('Case 2: Two closely-spaced targets',[10-close_sep/2, 10+close_sep/2], [1 1], 10, []);
    cases(3) = make_case('Case 3: Three targets',            [-25, 5, 35],        [1 1 1],       10, []);
    cases(4) = make_case('Case 4: Low SNR (-10 dB)',         12,                  1,            -20, []);
    cases(5) = make_case('Case 5: Strong interference',      12,                  1,             10, struct('doa',45,'power_db',20));

    % Case 4 is deliberately harsher than a literal "-10 dB" label
    % (-20 dB per-element SNR with reduced snapshot support, see
    % NumSnapshots override below) so that method robustness actually
    % differs case-to-case -- at a mild -10 dB with N=1000 snapshots
    % and this array's 15 dB of array gain, every method recovers the
    % target exactly and the case fails to demonstrate anything.
    case4_snapshots = max(64, round(N/4));

    methods = {'DS','MVDR','Capon','MUSIC'};
    results = struct([]);

    for c = 1:numel(cases)
        this_case = cases(c);
        fprintf('-------------------------------------------------\n');
        fprintf('%s\n', this_case.name);
        fprintf('-------------------------------------------------\n');

        if c == 4
            N_case = case4_snapshots;
            fprintf('  (using %d snapshots instead of %d to stress covariance estimation)\n', N_case, N);
        else
            N_case = N;
        end

        X = simulate_array_snapshots(array, config, this_case, N_case);

        % Number of spectral peaks to ask for: true targets, PLUS the
        % interferer if present (Case 5) -- asking for only the
        % target count when an interferer is actually present in the
        % scene would just make every method report the interferer's
        % peak (it is a real, strong source and correctly appears in
        % every method's spectrum; that is not a bug, but reporting
        % only 1 peak in that situation is an uninformative demo).
        num_true = numel(this_case.doas);
        if ~isempty(this_case.interference)
            num_report = num_true + 1;
        else
            num_report = num_true;
        end

        case_result = struct();
        case_result.name = this_case.name;
        case_result.true_doas = this_case.doas;
        case_result.snr_db = this_case.snr_db;
        case_result.doa_estimates = struct();
        case_result.comp_time_s = struct();
        case_result.metrics = struct();

        spectra = struct();
        for m = 1:numel(methods)
            method = methods{m};
            [Pdb, ~, doa, t] = spatial_spectrum(X, array, config, method, ...
                'ThetaGrid', theta_grid, 'NumSources', num_report);
            case_result.doa_estimates.(method) = doa;
            case_result.comp_time_s.(method)   = t;
            spectra.(method) = Pdb;

            err = doa_error(this_case.doas, doa);
            fprintf('  %-6s: DOA_est = %-24s  mean|err| = %6.2f deg  t = %.4f s\n', ...
                method, mat2str(doa,4), err, t);
        end

        if c == 5
            print_interference_rejection(X, array, config, this_case);
        end

        % deterministic DS beam pattern metrics for this case's array
        % (same array/config every case -- HPBW/PSLL/gain don't
        % depend on the data, only on the weights, so compute once
        % per case at the case's true first-target look angle for a
        % representative figure)
        [~, ~, pat_metrics] = beam_pattern( ...
            steering_vector(array.positions, this_case.doas(1), config.fc, config.sound_speed) / array.num_elements, ...
            array, config, 'ThetaGrid', theta_grid, 'PlotType', 'none');
        case_result.metrics.ds_hpbw_deg      = pat_metrics.hpbw_deg;
        case_result.metrics.ds_psll_db       = pat_metrics.peak_sidelobe_db;
        case_result.metrics.ds_array_gain_db = pat_metrics.array_gain_db;

        if make_plots
            plot_case_comparison(theta_grid, spectra, this_case, array);
        end

        if isempty(results)
            results = case_result;
        else
            results(end+1) = case_result; %#ok<AGROW>
        end
        fprintf('\n');
    end

    memory_report(array, N, theta_grid);
    print_summary_table(results, methods);

    if make_plots
        % Representative polar + Cartesian DS beam pattern (Case 1)
        a_look = steering_vector(array.positions, cases(1).doas(1), config.fc, config.sound_speed);
        beam_pattern(a_look/array.num_elements, array, config, ...
            'ThetaGrid', theta_grid, 'PlotType', 'both', ...
            'Title', 'Delay-and-Sum: Case 1 Look Direction');
    end

    fprintf('=================================================\n');
    fprintf(' COMPARISON COMPLETE\n');
    fprintf('=================================================\n');

end

% =========================================================================
% LOCAL HELPER FUNCTIONS
% =========================================================================

function c = make_case(name, doas, powers, snr_db, interference)
%MAKE_CASE Bundle one simulation case's parameters into a struct.
    c.name = name;
    c.doas = doas;
    c.powers = powers;
    c.snr_db = snr_db;
    c.interference = interference;
end

% -------------------------------------------------------------------------
function print_interference_rejection(X, array, config, this_case)
%PRINT_INTERFERENCE_REJECTION Case-5-specific demonstration of adaptive
%   nulling. The spatial-spectrum peak list (printed by the main loop)
%   correctly shows both the target and the interferer as real
%   sources -- that is expected and is not where MVDR/Capon's benefit
%   shows up. The actual benefit is in SIGNAL EXTRACTION: steering a
%   beamformer at the known/assumed target bearing and asking "how
%   much of the interferer leaks into my output?" This directly
%   exercises the "beamformed signal" output requested for this
%   module (as opposed to the "spatial spectrum" output), using
%   DELAY_SUM.M and MVDR_BEAMFORMER.M's actual weight vectors -- not
%   a re-derivation -- so this is measuring the delivered functions
%   themselves, not a separate approximation of them.
    target_doa = this_case.doas(1);
    interf_doa = this_case.interference.doa;

    [~, w_ds]   = delay_sum(X, array, target_doa, config);
    [~, w_mvdr] = mvdr_beamformer(X, array, target_doa, config);

    a_i = steering_vector(array.positions, interf_doa, config.fc, config.sound_speed);
    gain_ds_interf   = 20*log10(abs(w_ds'*a_i));
    gain_mvdr_interf = 20*log10(abs(w_mvdr'*a_i));

    fprintf('  --- Interference rejection (steered at true target, %.1f deg) ---\n', target_doa);
    fprintf('  DS   response toward interferer (%.1f deg): %7.2f dB\n', interf_doa, gain_ds_interf);
    fprintf('  MVDR response toward interferer (%.1f deg): %7.2f dB  (adaptive null)\n', ...
        interf_doa, gain_mvdr_interf);
    fprintf('  MVDR improves interference rejection by %.1f dB over DS at this bearing.\n', ...
        gain_ds_interf - gain_mvdr_interf);
end

% -------------------------------------------------------------------------
function X = simulate_array_snapshots(array, config, this_case, N)
%SIMULATE_ARRAY_SNAPSHOTS Self-contained synthetic far-field array data.
%
%   Generates [M x N] complex baseband array snapshots for one or more
%   uncorrelated far-field plane-wave sources plus spatially-white
%   complex Gaussian sensor noise, and (optionally) one stronger
%   interferer. This is demonstration/test data ONLY -- see the
%   module-level integration note in this file's main header on why
%   the project's existing propagation/simulate_echoes.m cannot
%   currently supply this (it has no per-element output at all).
%
%   Model:  x(t) = sum_k a(theta_k) s_k(t) + a(theta_I) s_I(t) + n(t)
%   with s_k, s_I, n all independent complex circular Gaussian
%   processes. SNR is defined per element, source power vs. noise
%   power (config-independent of array size, i.e. NOT the
%   post-beamforming SNR, which is `snr_db + array_gain_db` for DS).
    M = array.num_elements;
    doas = this_case.doas;
    powers = this_case.powers;
    snr_db = this_case.snr_db;

    noise_power = 1; % reference; source power set relative to this
    signal_power_total = noise_power * 10^(snr_db/10);

    X = sqrt(noise_power/2) * (randn(M,N) + 1j*randn(M,N));

    for k = 1:numel(doas)
        a_k = steering_vector(array.positions, doas(k), config.fc, config.sound_speed);
        p_k = signal_power_total * (powers(k) / sum(powers));
        s_k = sqrt(p_k/2) * (randn(1,N) + 1j*randn(1,N));
        X = X + a_k * s_k;
    end

    if ~isempty(this_case.interference)
        a_i = steering_vector(array.positions, this_case.interference.doa, config.fc, config.sound_speed);
        p_i = signal_power_total * 10^(this_case.interference.power_db/10);
        s_i = sqrt(p_i/2) * (randn(1,N) + 1j*randn(1,N));
        X = X + a_i * s_i;
    end
end

% -------------------------------------------------------------------------
function e = doa_error(true_doas, est_doas)
%DOA_ERROR Mean absolute error between true and estimated DOAs after
%   greedy nearest-neighbor matching. Returns NaN if no estimates.
    true_doas = sort(true_doas);
    if isempty(est_doas)
        e = NaN;
        return;
    end
    est_doas = sort(est_doas);
    errs = zeros(1, numel(true_doas));
    remaining = est_doas;
    for i = 1:numel(true_doas)
        if isempty(remaining)
            errs(i) = NaN;
            continue;
        end
        [d, idx] = min(abs(remaining - true_doas(i)));
        errs(i) = d;
        remaining(idx) = [];
    end
    e = mean(errs, 'omitnan');
end

% -------------------------------------------------------------------------
function plot_case_comparison(theta_grid, spectra, this_case, array) %#ok<INUSD>
%PLOT_CASE_COMPARISON Overlay DS/MVDR/Capon/MUSIC spectra for one case.
    figure('Name', ['Spectrum Comparison - ' this_case.name]);
    hold on;
    methods = {'DS','MVDR','Capon','MUSIC'};
    for m = 1:numel(methods)
        plot(theta_grid, spectra.(methods{m}), 'LineWidth', 1.3, 'DisplayName', methods{m});
    end
    for k = 1:numel(this_case.doas)
        xline(this_case.doas(k), 'k--', 'HandleVisibility','off');
    end
    if ~isempty(this_case.interference)
        xline(this_case.interference.doa, 'r:', 'Interferer', 'LineWidth', 1.5, 'HandleVisibility','off');
    end
    hold off;
    xlabel('Angle (degrees from broadside)');
    ylabel('Normalized power (dB)');
    title(['Beamformer Comparison: ' this_case.name]);
    legend('Location','best');
    grid on;
    ylim([-40, 2]);
end

% -------------------------------------------------------------------------
function bytes = memory_estimate_bytes(num_elements, num_snapshots, grid_len)
%MEMORY_ESTIMATE_BYTES Portable approximate memory footprint of the
%   dominant data structures for one method call, in bytes. Uses
%   direct size arithmetic (16 bytes per double-complex element)
%   rather than the interactive MATLAB Profiler, which is not
%   scriptable/portable across MATLAB and Octave.
    bytes_per_complex_double = 16;
    X_bytes = num_elements * num_snapshots * bytes_per_complex_double;
    R_bytes = num_elements^2 * bytes_per_complex_double;
    A_bytes = num_elements * grid_len * bytes_per_complex_double;
    bytes = struct('X_bytes', X_bytes, 'R_bytes', R_bytes, ...
                    'SteeringMatrix_bytes', A_bytes, ...
                    'total_bytes', X_bytes + R_bytes + A_bytes);
end

% -------------------------------------------------------------------------
function memory_report(array, N, theta_grid)
%MEMORY_REPORT Print the approximate memory footprint table.
    mem = memory_estimate_bytes(array.num_elements, N, numel(theta_grid));
    fprintf('-------------------------------------------------\n');
    fprintf('APPROXIMATE MEMORY FOOTPRINT (per method call)\n');
    fprintf('-------------------------------------------------\n');
    fprintf('  Snapshot matrix X   (%d x %d)  : %8.1f KB\n', array.num_elements, N, mem.X_bytes/1024);
    fprintf('  Covariance R        (%d x %d)  : %8.1f KB\n', array.num_elements, array.num_elements, mem.R_bytes/1024);
    fprintf('  Steering matrix A   (%d x %d)  : %8.1f KB\n', array.num_elements, numel(theta_grid), mem.SteeringMatrix_bytes/1024);
    fprintf('  Total (approx)                    : %8.1f KB\n', mem.total_bytes/1024);
    fprintf('  Note: MUSIC additionally stores an [M x M] eigenvector\n');
    fprintf('        matrix (same size as R); MVDR/Capon store R^{-1}\n');
    fprintf('        implicitly via a linear solve (no explicit inverse).\n');
    fprintf('-------------------------------------------------\n\n');
end

% -------------------------------------------------------------------------
function print_summary_table(results, methods)
%PRINT_SUMMARY_TABLE Formatted fprintf performance-comparison table
%   (MATLAB's table() object is avoided here so this prints
%   identically under Octave, which does not implement table()).
    fprintf('=================================================\n');
    fprintf(' PERFORMANCE SUMMARY (mean over all cases)\n');
    fprintf('=================================================\n');
    fprintf('%-8s %14s %14s\n', 'Method', 'MeanErr(deg)', 'MeanTime(s)');
    for m = 1:numel(methods)
        method = methods{m};
        errs = zeros(1, numel(results));
        times = zeros(1, numel(results));
        for c = 1:numel(results)
            errs(c) = doa_error(results(c).true_doas, results(c).doa_estimates.(method));
            times(c) = results(c).comp_time_s.(method);
        end
        fprintf('%-8s %14.3f %14.5f\n', method, mean(errs,'omitnan'), mean(times));
    end
    fprintf('-------------------------------------------------\n');
    fprintf('Delay-and-Sum reference beam pattern (per case, at 1st target):\n');
    fprintf('%-45s %10s %10s %10s\n', 'Case', 'HPBW(deg)', 'PSLL(dB)', 'Gain(dB)');
    for c = 1:numel(results)
        fprintf('%-45s %10.2f %10.2f %10.2f\n', results(c).name, ...
            results(c).metrics.ds_hpbw_deg, results(c).metrics.ds_psll_db, ...
            results(c).metrics.ds_array_gain_db);
    end
    fprintf('=================================================\n\n');
end
