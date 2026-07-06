function sp = sound_speed_profile(water_props, z_query, useSpline)
%SOUND_SPEED_PROFILE Interpolate c(z) and provide effective constant sound speed.
%
%   sp = sound_speed_profile(water_props, z_query)
%   sp = sound_speed_profile(water_props, z_query, useSpline)
%
% Inputs:
%   water_props - struct from water_properties.m (fields: z, c).
%   z_query     - query depths (m), vector.
%   useSpline   - optional logical (default true). If false, piecewise-linear.
%
% Output struct:
%   sp.c_local - interpolated local sound speed at z_query.
%   sp.c_eff   - harmonic-path-averaged effective sound speed.

validateattributes(water_props, {'struct'},{'nonempty'}, mfilename,'water_props',1);
validateattributes(z_query, {'numeric'},{'real','finite','vector','nonempty'}, mfilename,'z_query',2);
if nargin < 3
    useSpline = true;
end
validateattributes(useSpline, {'numeric','logical'},{'scalar','finite'}, mfilename,'useSpline',3);
useSpline = logical(useSpline);

req = {'z','c'};
for k = 1:numel(req)
    if ~isfield(water_props, req{k})
        error('%s: water_props missing field "%s".', mfilename, req{k});
    end
end

z_tab = water_props.z(:);
c_tab = water_props.c(:);
validateattributes(z_tab, {'numeric'},{'real','finite','vector','nonempty'}, mfilename,'water_props.z',1);
validateattributes(c_tab, {'numeric'},{'real','finite','vector','nonempty'}, mfilename,'water_props.c',1);
if numel(z_tab) ~= numel(c_tab)
    error('%s: water_props.z and water_props.c must have the same length.', mfilename);
end
if any(diff(z_tab) <= 0)
    error('%s: water_props.z must be strictly monotonically increasing.', mfilename);
end

zq = z_query(:);

% Interpolate; allow extrapolation as defensively filled.
if useSpline && numel(z_tab) >= 4
    c_local = interp1(z_tab, c_tab, zq, 'spline', 'extrap');
else
    c_local = interp1(z_tab, c_tab, zq, 'linear', 'extrap');
end

% Effective harmonic mean over tabulated profile.
% sp.c_eff = 1 / mean(1./c_tab)
if any(c_tab <= 0)
    error('%s: sound speed must be positive everywhere.', mfilename);
end
c_eff = 1 ./ mean(1 ./ c_tab);

sp.c_local = reshape(c_local, size(z_query));
sp.c_eff   = c_eff;

end

