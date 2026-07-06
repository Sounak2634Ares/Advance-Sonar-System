function out = transmission_loss(range_m_vec, alpha_dB_per_km, c)
%TRANSMISSION_LOSS Combine geometric spreading and absorption with range clamping.
%
%   out = transmission_loss(range_m_vec, alpha_dB_per_km, c)
%
% Inputs:
%   range_m_vec      - Range vector/matrix (meters), any shape.
%   alpha_dB_per_km  - Absorption coefficient at fc (dB/km), scalar.
%   c                - Sound speed (m/s), scalar (defensive; not used in TL itself).
%
% Output struct:
%   out.TL_dB  - Two-way transmission loss in dB, same size as range_m_vec.
%   out.G_tvg  - Two-way TVG gain (linear amplitude) = 10^(TL/20), same size.

validateattributes(range_m_vec, {'numeric'},{'real','finite','nonempty'}, mfilename,'range_m_vec',1);
validateattributes(alpha_dB_per_km, {'numeric'},{'real','finite','scalar','nonnegative'}, mfilename,'alpha_dB_per_km',2);
validateattributes(c, {'numeric'},{'real','finite','scalar','positive'}, mfilename,'c',3); %#ok<NASGU>

eps0 = eps('double');
R = range_m_vec;

% CRITICAL clamp to prevent log10(0) = -Inf
R_safe = max(R, eps0);

TL_spreading = 40 .* log10(R_safe);
TL_absorption = 2 .* alpha_dB_per_km .* (R_safe./1000);

TL_2way = TL_spreading + TL_absorption;

G_tvg = 10.^(TL_2way./20);

out.TL_dB = TL_2way;
out.G_tvg = G_tvg;

end

