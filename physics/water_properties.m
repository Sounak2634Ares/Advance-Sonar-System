function props = water_properties(T_profile, S_profile, z_profile, phi_lat)
%WATER_PROPERTIES Convert oceanographic profiles to sound speed and density fields.
%
%   props = water_properties(T_profile, S_profile, z_profile, phi_lat)
%
% Inputs:
%   T_profile  - Temperature profile (°C). Vector, same length as z_profile.
%   S_profile  - Salinity profile (ppt or PSU). Vector, same length as z_profile.
%   z_profile  - Depth profile (m). Vector, same length as T_profile.
%   phi_lat    - Latitude for gravity correction (deg). Optional; currently unused.
%
% Output (struct):
%   props.c   - Sound speed profile (m/s), same shape as inputs (column).
%   props.rho - Linearized density profile (kg/m^3), same shape as props.c.
%   props.T, props.S, props.z - sanitized input vectors (column).

validateattributes(T_profile, {'numeric'},{'real','finite','vector','nonempty'}, mfilename,'T_profile',1);
validateattributes(S_profile, {'numeric'},{'real','finite','vector','nonempty'}, mfilename,'S_profile',2);
validateattributes(z_profile, {'numeric'},{'real','finite','vector','nonempty'}, mfilename,'z_profile',3);
if nargin < 4
    phi_lat = NaN;
end
validateattributes(phi_lat, {'numeric'},{'real','finite','scalar'}, mfilename,'phi_lat',4);

T = T_profile(:);
S = S_profile(:);
z = z_profile(:);

n = numel(z);
if numel(T) ~= n || numel(S) ~= n
    error('%s: T_profile, S_profile, z_profile must have the same length.', mfilename);
end

% Monotonic increasing depth requirement
if any(diff(z) <= 0)
    error('%s: z_profile must be strictly monotonically increasing.', mfilename);
end

% Mackenzie sound speed equation (as specified)
props.c = 1448.96 + 4.591*T - 5.304e-2*T.^2 + 2.374e-4*T.^3 + 1.340*(S-35) ...
          + 1.630e-2*z + 1.675e-7*z.^2 - 1.025e-2*T.*(S-35) - 7.139e-13*T.*z.^3;

% Density: linearized UNESCO/EOS-80 simplified (as specified)
% rho0[kg/m^3] and beta_T, beta_S not explicitly given in prompt; use typical
% small-signal coefficients consistent with linearized form.
% If you later supply beta_T/beta_S, replace these values.
rho0 = 1027;
T0 = 0;   % reference temperature (°C) - consistent with T-T0 term
S0 = 35;  % reference salinity (ppt)

% Typical compressibility/thermal expansion coefficients (approx)
% so that rho decreases with increasing T and increases with increasing S.
beta_T = 2.6e-4;
beta_S = 7.6e-4;

props.rho = rho0 .* (1 - beta_T*(T - T0) + beta_S*(S - S0));

props.T = T;
props.S = S;
props.z = z;

end

