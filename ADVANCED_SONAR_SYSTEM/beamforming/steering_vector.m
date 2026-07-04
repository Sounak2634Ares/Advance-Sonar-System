function A = steering_vector(positions, theta_deg, frequency, sound_speed, phi_deg)
%STEERING_VECTOR Far-field narrowband array manifold vector(s).
%
% -------------------------------------------------------------------
% PURPOSE
%   Computes the complex array manifold (steering) vector(s) that
%   every beamformer in this module (DELAY_SUM, MVDR_BEAMFORMER,
%   CAPON_BEAMFORMER, MUSIC_DOA, BEAM_PATTERN, SPATIAL_SPECTRUM) is
%   built on. Implemented once, generically, in 3-D so that Circular
%   and Planar arrays can be added later purely by changing
%   ARRAY_GEOMETRY.M -- this file never needs to change.
%
% MATHEMATICAL THEORY
%   For a far-field plane wave arriving from azimuth theta (measured
%   from broadside) and elevation phi (measured from the array's
%   x-y plane), the unit vector pointing toward the source is
%
%       k_hat(theta,phi) = [ sin(theta)cos(phi) ;
%                             cos(theta)cos(phi) ;
%                             sin(phi)          ]
%
%   The phase of the wavefront at an element located at p_n (meters,
%   relative to the array phase center) relative to the phase center
%   is  -2*pi/lambda * (p_n . k_hat), so the steering vector is
%
%       a_n(theta,phi) = exp( -j * 2*pi/lambda * (p_n . k_hat) )   (1)
%
%   with lambda = sound_speed / frequency (narrowband assumption --
%   see Known Limitations below).
%
%   Equivalence to time-delay ("delay-and-sum"):
%       The geometric propagation delay at element n relative to the
%       phase center is  tau_n = -(p_n . k_hat)/c.  For a NARROWBAND
%       signal at frequency f, delaying by tau_n is exactly equivalent
%       to multiplying by exp(-j*2*pi*f*tau_n), which is precisely
%       Eq.(1). This is why phase-shift beamforming and "delay-and-
%       sum" beamforming are the same operation for narrowband arrays
%       (Van Trees, Optimum Array Processing, Sec. 2.3).
%
%   For a Uniform Linear Array with elements on the local x-axis
%   (phi = 0, the standard case used by ARRAY_GEOMETRY.M's 'ULA'
%   type), Eq.(1) reduces to the familiar 1-D form
%
%       a_n(theta) = exp( -j * 2*pi/lambda * d_n * sin(theta) )    (2)
%
%   where d_n is the element's position along the array line.
%
%   IMPORTANT INTEGRATION NOTE: the pre-existing
%   beamforming/multibeam_configuration.m in this project computes a
%   steering_matrix using "exp(-1j*2*pi*element_positions'*sin(theta))"
%   -- i.e. Eq.(2) WITHOUT the 1/lambda factor. That is only
%   dimensionally correct if lambda = 1 m, which it is not for this
%   system (lambda = 0.03 m at fc = 50 kHz). This function implements
%   the physically correct Eq.(1)/(2) and intentionally does NOT reuse
%   multibeam_configuration.m's steering_matrix. See the module-level
%   integration notes for a recommendation.
%
% ALGORITHM
%   1. Build the [K x 3] matrix of unit direction vectors k_hat for
%      the K requested (theta,phi) pairs.
%   2. A = exp( -1j * (2*pi/lambda) * positions * k_hat' )   -- a
%      single [M x 3] * [3 x K] matrix multiply gives the full
%      [M x K] steering matrix for all K angles at once (vectorized;
%      no per-angle loop).
%
% INPUT PARAMETERS
%   positions   : [M x 3] element coordinates in meters (from
%                 ARRAY_GEOMETRY.M), column order (x,y,z).
%   theta_deg   : azimuth angle(s) in degrees, measured from
%                 broadside. Scalar or vector of length K.
%   frequency   : operating (narrowband/center) frequency in Hz,
%                 scalar. Typically config.fc.
%   sound_speed : propagation speed in m/s, scalar. Typically
%                 config.sound_speed.
%   phi_deg     : (optional) elevation angle(s) in degrees, default
%                 0. Must be scalar or the same length as theta_deg.
%                 A 1-D ULA cannot resolve elevation independently of
%                 azimuth (cone-of-ambiguity); this argument exists so
%                 future circular/planar geometries can use the same
%                 function unmodified.
%
% OUTPUT PARAMETERS
%   A : [M x K] complex steering matrix, one column per requested
%       angle, each column with unit-modulus entries (|A(:,k)|=1).
%
% COMPLEXITY
%   Time  : O(M*K) for the matrix multiply and complex exponential.
%   Memory: O(M*K) for the output matrix (dominant cost for large
%           angle-grid sweeps, e.g. M=32, K=721 => ~23k complex
%           doubles = ~370 KB -- negligible for any realistic array).
%
% ASSUMPTIONS
%   - Narrowband assumption: a single "frequency" is used for the
%     whole signal bandwidth. This project's chirp uses fc = 50 kHz
%     over a 60 kHz bandwidth (fractional bandwidth ~120%), which is
%     genuinely wideband; using fc as a single reference frequency is
%     the standard simplification for classical DS/MVDR/Capon/MUSIC
%     and is what is implemented here. See "Known Limitations" in the
%     module-level notes for a wideband (subband/frequency-domain)
%     extension path.
%   - Far-field (plane-wave) propagation.
%
% REFERENCES
%   [1] H. L. Van Trees, "Optimum Array Processing," Wiley, 2002.
%   [2] D. H. Johnson & D. E. Dudgeon, "Array Signal Processing,"
%       Prentice Hall, 1993.
%
% See also ARRAY_GEOMETRY, DELAY_SUM, MVDR_BEAMFORMER
% -------------------------------------------------------------------

    if nargin < 5 || isempty(phi_deg)
        phi_deg = 0;
    end

    validateattributes(positions, {'numeric'}, {'2d','ncols',3,'real'}, ...
        'steering_vector', 'positions');
    validateattributes(theta_deg, {'numeric'}, {'real'}, ...
        'steering_vector', 'theta_deg');
    validateattributes(frequency, {'numeric'}, {'scalar','real','positive'}, ...
        'steering_vector', 'frequency');
    validateattributes(sound_speed, {'numeric'}, {'scalar','real','positive'}, ...
        'steering_vector', 'sound_speed');

    theta_deg = theta_deg(:)';   % 1 x K row
    K = numel(theta_deg);

    if isscalar(phi_deg)
        phi_deg = repmat(phi_deg, 1, K);
    else
        phi_deg = phi_deg(:)';
        if numel(phi_deg) ~= K
            error('steering_vector:SizeMismatch', ...
                'phi_deg must be scalar or the same length as theta_deg.');
        end
    end

    lambda = sound_speed / frequency;

    % [3 x K] unit direction vectors: rows = (x,y,z) components
    k_hat = [ sind(theta_deg) .* cosd(phi_deg) ; ...
              cosd(theta_deg) .* cosd(phi_deg) ; ...
              sind(phi_deg) ];

    % [M x 3] * [3 x K] = [M x K], single vectorized matrix multiply
    A = exp(-1j * (2*pi/lambda) * (positions * k_hat));

end
