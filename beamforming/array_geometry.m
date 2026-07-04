function array = array_geometry(config, varargin)
%ARRAY_GEOMETRY Generate sonar array element positions and metadata.
%
% -------------------------------------------------------------------
% PURPOSE
%   Produces a self-contained array description (element positions,
%   spacing, aperture) consumed by STEERING_VECTOR.M and all four
%   beamformer modules (DELAY_SUM, MVDR_BEAMFORMER, CAPON_BEAMFORMER,
%   MUSIC_DOA). This is the single point of extension for new array
%   topologies: adding a circular or planar array means adding a new
%   CASE below and nowhere else -- STEERING_VECTOR.M is fully generic
%   3-D geometry and never needs to change.
%
% MATHEMATICAL THEORY
%   Uniform Linear Array (ULA):
%       Elements lie on a line (taken here as the local x-axis) at
%           p_n = n * d,   n = 0, 1, ..., M-1
%       and are then centered so the array's phase center is the
%       origin:
%           p_n <- p_n - mean(p)
%       Each element position is stored as a 3-D coordinate (x,y,z)
%       with y = z = 0 for a ULA, which is what allows
%       STEERING_VECTOR.M to use one general dot-product formula for
%       any future geometry (circular, planar) without modification.
%
%   Design rule-of-thumb (grating lobes):
%       To avoid spatial aliasing (grating lobes) over the full
%       visible region, spacing should satisfy d <= lambda/2 at the
%       highest frequency of interest. This function warns if that
%       is violated for config.fc.
%
% INPUT PARAMETERS
%   config      : struct, must contain (unless overridden, see below):
%                   .array_elements   number of hydrophones (M)
%                   .element_spacing  inter-element spacing, meters
%                   .fc               operating frequency, Hz  (for
%                                     the grating-lobe sanity check)
%                   .sound_speed      m/s (for the same check)
%
%   Name-Value overrides (optional, useful for unit tests / sweeps
%   without needing a full config struct):
%       'ArrayType'   : 'ULA' (default). 'circular' and 'planar' are
%                       recognized names reserved for future work and
%                       currently raise a clear "not yet implemented"
%                       error rather than silently returning wrong
%                       geometry.
%       'NumElements' : overrides config.array_elements
%       'Spacing'     : overrides config.element_spacing
%
% OUTPUT PARAMETERS
%   array : struct with fields
%       .type              'ULA'
%       .num_elements       M
%       .spacing            d (m)
%       .positions          [M x 3] element coordinates (m), centered
%       .aperture           physical length of the array (m)
%       .config_fc          fc used for the grating-lobe check (Hz)
%
% COMPLEXITY
%   Time  : O(M)
%   Memory: O(M) -- a [M x 3] position matrix, negligible for any
%           realistic array size.
%
% ASSUMPTIONS
%   - ULA elements are collinear and equally spaced (no per-element
%     position error / manufacturing tolerance modeling).
%   - The array's local x-axis is the array line; broadside is the
%     local y-axis. How this local frame maps onto the vehicle/world
%     (x,y,z) frame used elsewhere in the project (see
%     targets/generate_targets.m) is not yet defined anywhere in the
%     existing codebase (no boresight/orientation parameter exists in
%     config). That reconciliation is out of scope for this module;
%     see the integration notes delivered alongside this file.
%
% REFERENCES
%   [1] H. L. Van Trees, "Optimum Array Processing," Wiley, 2002, Ch.2.
%   [2] D. H. Johnson & D. E. Dudgeon, "Array Signal Processing,"
%       Prentice Hall, 1993, Ch.2-3.
%
% See also STEERING_VECTOR, DELAY_SUM, MVDR_BEAMFORMER
% -------------------------------------------------------------------

    p = inputParser;
    addRequired(p, 'config', @(c) isstruct(c) || isempty(c));
    addParameter(p, 'ArrayType',   'ULA', @(s) ischar(s) || isstring(s));
    addParameter(p, 'NumElements', [],    @(v) isempty(v) || (isnumeric(v) && isscalar(v)));
    addParameter(p, 'Spacing',     [],    @(v) isempty(v) || (isnumeric(v) && isscalar(v)));
    parse(p, config, varargin{:});
    opts = p.Results;

    array_type = lower(char(opts.ArrayType));

    % ---------------------------------------------------------------
    % Resolve parameters: explicit Name-Value overrides win, else
    % pull from config, matching the project-wide "config struct"
    % convention (see config/system_config.m).
    % ---------------------------------------------------------------
    if ~isempty(opts.NumElements)
        M = opts.NumElements;
    else
        if ~isstruct(config) || ~isfield(config, 'array_elements')
            error('array_geometry:MissingParameter', ...
                ['config.array_elements not found and no NumElements ' ...
                 'override supplied.']);
        end
        M = config.array_elements;
    end

    if ~isempty(opts.Spacing)
        d = opts.Spacing;
    else
        if ~isstruct(config) || ~isfield(config, 'element_spacing')
            error('array_geometry:MissingParameter', ...
                ['config.element_spacing not found and no Spacing ' ...
                 'override supplied.']);
        end
        d = config.element_spacing;
    end

    % ---------------------------------------------------------------
    % Validate
    % ---------------------------------------------------------------
    validateattributes(M, {'numeric'}, {'scalar','integer','positive'}, ...
        'array_geometry', 'NumElements');
    validateattributes(d, {'numeric'}, {'scalar','real','positive'}, ...
        'array_geometry', 'Spacing');

    switch array_type
        case 'ula'
            positions = (0:M-1)' * d;              % [M x 1], local x-axis
            positions = positions - mean(positions); % center phase reference
            positions3 = [positions, zeros(M,1), zeros(M,1)]; % [M x 3]
            aperture = (M-1) * d;

        case 'circular'
            error('array_geometry:NotImplemented', ...
                ['Circular array geometry is reserved for a future ' ...
                 'phase and is not yet implemented. STEERING_VECTOR.M ' ...
                 'already supports arbitrary 3-D positions, so adding ' ...
                 'this case is the only change required later.']);

        case 'planar'
            error('array_geometry:NotImplemented', ...
                ['Planar array geometry is reserved for a future ' ...
                 'phase and is not yet implemented. STEERING_VECTOR.M ' ...
                 'already supports arbitrary 3-D positions, so adding ' ...
                 'this case is the only change required later.']);

        otherwise
            error('array_geometry:UnknownType', ...
                'Unknown ArrayType "%s". Supported now: ULA.', array_type);
    end

    % ---------------------------------------------------------------
    % Grating-lobe design check (warning only, never fatal -- the
    % array may intentionally be sparse for a lower operating band)
    % ---------------------------------------------------------------
    fc = NaN; c = NaN;
    if isstruct(config)
        if isfield(config, 'fc'),          fc = config.fc;          end
        if isfield(config, 'sound_speed'), c  = config.sound_speed; end
    end
    if isfinite(fc) && isfinite(c)
        lambda = c / fc;
        if d > lambda/2 + 1e-12
            warning('array_geometry:GratingLobeRisk', ...
                ['Element spacing (%.5f m) exceeds lambda/2 (%.5f m) at ' ...
                 'fc = %.0f Hz. Grating lobes may appear in the visible ' ...
                 'region for steering angles away from broadside.'], ...
                 d, lambda/2, fc);
        end
    end

    % ---------------------------------------------------------------
    % Assemble output
    % ---------------------------------------------------------------
    array = struct();
    array.type          = upper(array_type);
    array.num_elements   = M;
    array.spacing        = d;
    array.positions       = positions3;
    array.aperture        = aperture;
    array.config_fc       = fc;

    fprintf('-----------------------------------\n');
    fprintf('ARRAY GEOMETRY (%s)\n', array.type);
    fprintf('-----------------------------------\n');
    fprintf('Elements        : %d\n', M);
    fprintf('Spacing         : %.5f m\n', d);
    fprintf('Aperture        : %.4f m\n', aperture);
    if isfinite(fc) && isfinite(c)
        fprintf('d / lambda      : %.3f (at fc = %.0f Hz)\n', d/(c/fc), fc);
    end
    fprintf('-----------------------------------\n');

end
