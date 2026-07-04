function [pattern_db, theta_grid, metrics] = beam_pattern(w, array, config, varargin)
%BEAM_PATTERN Deterministic array beam pattern, metrics, and plots.
%
% -------------------------------------------------------------------
% PURPOSE
%   Computes the radiation (response) pattern |w^H a(theta)|^2 for a
%   FIXED weight vector w, swept over an angle grid, plus the three
%   deterministic performance figures the project asks for
%   (beamforming gain, resolution/HPBW, peak sidelobe level), and
%   optionally renders polar + Cartesian plots.
%
%   This complements SPATIAL_SPECTRUM.M: BEAM_PATTERN is for a weight
%   vector that is already fixed (e.g. DELAY_SUM's conventional
%   taper, or MVDR_BEAMFORMER's adapted weights frozen for one
%   snapshot), while SPATIAL_SPECTRUM/CAPON_BEAMFORMER/MUSIC_DOA
%   recompute their statistic fresh at every scanned angle from data.
%
% MATHEMATICAL THEORY
%   Normalized power pattern:
%       p(theta) = | w^H a(theta) |^2 / | w^H a(theta_0) |^2        (1)
%   where theta_0 is the look direction the weights were designed
%   for (taken here as the angle of peak response, so (1) is always
%   normalized to a 0 dB peak regardless of how w was produced).
%
%   Half-power (-3 dB) beamwidth (HPBW), i.e. "resolution": the
%   angular width of the mainlobe at the -3 dB level, found
%   numerically from p(theta) rather than only quoting the classical
%   broadside approximation
%       HPBW_broadside (rad) ~= 0.886 * lambda / (M*d)              (2)
%   (Eq. 2 is only exact at broadside with uniform, untapered
%   weights; the numerical HPBW in `metrics` is exact for whatever
%   w was actually supplied.)
%
%   Peak sidelobe level (PSLL): the height, in dB, of the highest
%   local maximum of p(theta) outside the mainlobe region (defined
%   here as +/- 1.5x the numerically-found HPBW around the peak).
%
%   Array (beamforming) gain: for spatially-white unit-power sensor
%   noise, the processing gain of weight vector w is
%       G = |w^H a(theta_0)|^2 / ( ||w||^2 )                        (3)
%   which reduces to G = M for the conventional w = a/M (matching
%   DELAY_SUM.M's `info.array_gain_db`), and is generally lower for
%   adapted (e.g. MVDR) weights, which trade some white-noise gain
%   for interference rejection.
%
% ALGORITHM
%   1. Build the steering matrix over theta_grid (one call).
%   2. p(theta) = |w^H A|^2, normalized to the peak (Eq. 1).
%   3. Find the peak, then numerically find the -3 dB crossings on
%      each side to get HPBW; exclude that region and find the
%      tallest remaining local max for PSLL.
%   4. Compute G via Eq. (3).
%   5. If plotting requested: render a polar and/or Cartesian pattern.
%
% INPUT PARAMETERS
%   w      : [M x 1] weight vector (a single column; for multiple
%            look angles, call once per column of a [M x K] weight
%            matrix -- kept as a single-pattern-per-call function so
%            plotting semantics stay simple and unambiguous).
%   array  : struct from ARRAY_GEOMETRY.M.
%   config : struct, needs .fc, .sound_speed.
%
%   Name-Value options:
%       'ThetaGrid' : [1 x K] angle grid, degrees. Default -90:0.1:90.
%       'PlotType'  : 'none' (default), 'polar', 'cartesian', or
%                     'both'.
%       'Title'     : string used in plot titles, default 'Beam Pattern'.
%
% OUTPUT PARAMETERS
%   pattern_db : [1 x K] normalized pattern in dB (0 dB peak).
%   theta_grid : the angle grid used (echoed back).
%   metrics    : struct with
%                 .hpbw_deg            half-power beamwidth, degrees
%                 .peak_sidelobe_db    PSLL relative to mainlobe peak
%                 .array_gain_db       10*log10(G), Eq.(3)
%                 .peak_angle_deg      angle of the pattern's peak
%
% COMPLEXITY
%   Time  : O(M*K) to build the steering matrix + O(K) for the
%           pattern and metric extraction.
%   Memory: O(M*K), dominated by the steering matrix for fine grids.
%
% ASSUMPTIONS
%   - Single weight vector per call (see INPUT PARAMETERS above).
%   - HPBW/PSLL extraction assumes a single dominant mainlobe; for
%     highly irregular adapted patterns (e.g. MVDR with multiple deep
%     nulls) the "mainlobe exclusion region" heuristic (+/-1.5x HPBW)
%     is a reasonable default but can occasionally include part of a
%     nearby lobe -- inspect the plot for adapted-weight patterns
%     rather than relying on PSLL alone.
%
% REFERENCES
%   [1] H. L. Van Trees, "Optimum Array Processing," Wiley, 2002, Ch.2.
%   [2] R. J. Mailloux, "Phased Array Antenna Handbook," Artech
%       House, 2nd ed., 2005 (beamwidth/sidelobe definitions).
%
% See also STEERING_VECTOR, SPATIAL_SPECTRUM, DELAY_SUM
% -------------------------------------------------------------------

    if nargin < 3
        error('beam_pattern:NotEnoughInputs', ...
            'Usage: [pattern_db,theta_grid,metrics] = beam_pattern(w, array, config, ...)');
    end

    p = inputParser;
    addParameter(p, 'ThetaGrid', -90:0.1:90, @(v) isnumeric(v) && isvector(v));
    addParameter(p, 'PlotType',  'none', @(s) any(strcmpi(s, {'none','polar','cartesian','both'})));
    addParameter(p, 'Title',     'Beam Pattern', @(s) ischar(s) || isstring(s));
    parse(p, varargin{:});
    theta_grid = p.Results.ThetaGrid(:)';
    plot_type  = lower(char(p.Results.PlotType));
    plot_title = char(p.Results.Title);

    w = w(:);
    if numel(w) ~= array.num_elements
        error('beam_pattern:DimensionMismatch', ...
            'Weight vector has %d entries but array has %d elements.', ...
            numel(w), array.num_elements);
    end

    A = steering_vector(array.positions, theta_grid, config.fc, config.sound_speed); % [M x K]
    resp = w' * A;                          % [1 x K]
    p_lin = abs(resp).^2;
    [peak_val, peak_idx] = max(p_lin);
    p_norm = p_lin / peak_val;
    pattern_db = 10*log10(max(p_norm, eps));

    % ---- HPBW: walk outward from the peak until dropping below -3 dB ----
    hpbw_deg = compute_hpbw(pattern_db, theta_grid, peak_idx);

    % ---- PSLL: exclude the mainlobe region, find tallest remaining lobe ----
    exclusion_halfwidth = 1.5 * max(hpbw_deg, mean(diff(theta_grid)));
    outside = abs(theta_grid - theta_grid(peak_idx)) > exclusion_halfwidth;
    if any(outside)
        peak_sidelobe_db = max(pattern_db(outside));
    else
        peak_sidelobe_db = -Inf; % grid too narrow to contain a sidelobe
    end

    % ---- array gain, Eq.(3) ----
    array_gain = peak_val / (w' * w);
    array_gain_db = 10*log10(real(array_gain));

    metrics = struct();
    metrics.hpbw_deg         = hpbw_deg;
    metrics.peak_sidelobe_db = peak_sidelobe_db;
    metrics.array_gain_db    = array_gain_db;
    metrics.peak_angle_deg   = theta_grid(peak_idx);

    if ~strcmp(plot_type, 'none')
        render_plots(pattern_db, theta_grid, metrics, plot_type, plot_title);
    end

end

% -----------------------------------------------------------------------
function hpbw_deg = compute_hpbw(pattern_db, theta_grid, peak_idx)
%COMPUTE_HPBW Numerically locate the -3 dB mainlobe width around the
%   peak by walking outward in both directions until the pattern
%   drops below -3 dB (linear interpolation between grid points for
%   sub-grid-step accuracy), or the grid edge is reached.
    K = numel(pattern_db);

    % walk left
    i = peak_idx;
    while i > 1 && pattern_db(i) >= -3
        i = i - 1;
    end
    if i == peak_idx
        left_deg = theta_grid(1);
    elseif pattern_db(i) < -3 && i < peak_idx
        left_deg = interp1(pattern_db(i:i+1), theta_grid(i:i+1), -3, 'linear');
    else
        left_deg = theta_grid(1);
    end

    % walk right
    j = peak_idx;
    while j < K && pattern_db(j) >= -3
        j = j + 1;
    end
    if j == peak_idx
        right_deg = theta_grid(end);
    elseif pattern_db(j) < -3 && j > peak_idx
        right_deg = interp1(pattern_db(j-1:j), theta_grid(j-1:j), -3, 'linear');
    else
        right_deg = theta_grid(end);
    end

    hpbw_deg = right_deg - left_deg;
end

% -----------------------------------------------------------------------
function render_plots(pattern_db, theta_grid, metrics, plot_type, plot_title)
%RENDER_PLOTS Polar and/or Cartesian beam pattern figures. Floors the
%   dB scale at -40 dB purely for display (does not affect metrics,
%   which are computed on the unfloored data beforehand).
    floor_db = -40;
    disp_db = max(pattern_db, floor_db);

    if any(strcmp(plot_type, {'cartesian','both'}))
        figure('Name', [plot_title ' (Cartesian)']);
        plot(theta_grid, disp_db, 'LineWidth', 1.5);
        hold on;
        yline(-3, '--', '-3 dB');
        xline(metrics.peak_angle_deg, ':', 'Peak');
        hold off;
        xlabel('Angle (degrees from broadside)');
        ylabel('Normalized power (dB)');
        title(sprintf('%s | HPBW=%.2f^\\circ, PSLL=%.1f dB, Gain=%.1f dB', ...
            plot_title, metrics.hpbw_deg, metrics.peak_sidelobe_db, metrics.array_gain_db));
        grid on;
        ylim([floor_db, 2]);
    end

    if any(strcmp(plot_type, {'polar','both'}))
        figure('Name', [plot_title ' (Polar)']);
        theta_rad = deg2rad(theta_grid);
        rho = disp_db - floor_db; % shift so floor maps to r=0
        polarplot(theta_rad, rho, 'LineWidth', 1.5);
        ax = gca;
        ax.ThetaZeroLocation = 'top';
        ax.ThetaDir = 'clockwise';
        title(sprintf('%s (Polar, floor at %d dB)', plot_title, floor_db));
    end
end
