function x_n_multipath = multipath_model(t, x_tx, x_horiz, zs, zr, H, rho1, c1, rho2, c2, targets_list, fc, alpha_dB_km, c)
%MULTIPATH_MODEL Multipath superposition using image-source method.
%
% This module extends the direct-path model by adding surface/bottom reflected arrivals.
%
% Inputs (defensive, flexible):
%   t             - time vector (s), [1 x Nt] or [Nt x 1]
%   x_tx         - transmit baseband signal samples, [1 x Nt] or [Nt x 1]
%   x_horiz      - horizontal range (m), scalar
%   zs, zr, H    - source depth, receiver depth, water depth (m)
%   rho1,c1       - densities/sound speed for medium 1 (water) (rho1 scalar, c1 scalar)
%   rho2,c2       - densities/sound speed for boundary medium 2 (bottom) (rho2 scalar, c2 scalar)
%   targets_list - struct array with fields:
%       * amp   (complex or real amplitude)
%       * range_multiplier (optional scalar to scale x_horiz for each target)
%       * phase0 (optional complex phase multiplier)
%   fc            - center frequency (Hz)
%   alpha_dB_km  - absorption coefficient (dB/km), scalar
%   c             - sound speed (m/s), scalar used for delays
%
% Output:
%   x_n_multipath - [Nt x numel(targets_list)] complex received signal per target.

validateattributes(t, {'numeric'},{'real','finite','vector','nonempty'}, mfilename,'t',1);
validateattributes(x_tx,{'numeric'},{'finite'}, mfilename,'x_tx',2);
if isempty(x_tx)
    x_n_multipath = zeros(numel(t(:)), numel(targets_list), 'like', complex(0));
    return;
end

validateattributes(x_horiz, {'numeric'},{'real','finite','scalar','nonnegative'}, mfilename,'x_horiz',3);
validateattributes(zs, {'numeric'},{'real','finite','scalar'}, mfilename,'zs',4);
validateattributes(zr, {'numeric'},{'real','finite','scalar'}, mfilename,'zr',5);
validateattributes(H, {'numeric'},{'real','finite','scalar','positive'}, mfilename,'H',6);
validateattributes(rho1, {'numeric'},{'real','finite','scalar','positive'}, mfilename,'rho1',7);
validateattributes(c1, {'numeric'},{'real','finite','scalar','positive'}, mfilename,'c1',8);
validateattributes(rho2, {'numeric'},{'real','finite','scalar','positive'}, mfilename,'rho2',9);
validateattributes(c2, {'numeric'},{'real','finite','scalar','positive'}, mfilename,'c2',10);
validateattributes(targets_list, {'struct'},{'nonempty'}, mfilename,'targets_list',11);
validateattributes(fc, {'numeric'},{'real','finite','scalar','positive'}, mfilename,'fc',12);
validateattributes(alpha_dB_km, {'numeric'},{'real','finite','scalar','nonnegative'}, mfilename,'alpha_dB_km',13);
validateattributes(c, {'numeric'},{'real','finite','scalar','positive'}, mfilename,'c',14);

% Ensure column time and signal
T = t(:);
Nt = numel(T);
if numel(x_tx) ~= Nt
    error('%s: length(x_tx) must equal length(t).', mfilename);
end
x_tx_col = x_tx(:);

M = numel(targets_list);
x_n_multipath = zeros(Nt, M, 'like', x_tx_col);

% Image path ranges (second-order set as specified)
% Use direct range and images for each target scaling.
for m = 1:M
    amp = 1;
    if isfield(targets_list(m),'amp')
        validateattributes(targets_list(m).amp, {'numeric'},{'finite'}, mfilename,'targets_list.amp',m);
        amp = targets_list(m).amp;
    end
    if isfield(targets_list(m),'phase0')
        amp = amp .* targets_list(m).phase0;
    end

    rm = 1;
    if isfield(targets_list(m),'range_multiplier')
        rm = targets_list(m).range_multiplier;
        validateattributes(rm, {'numeric'},{'real','finite','scalar','nonnegative'}, mfilename,'targets_list.range_multiplier',m);
    end
    xh = x_horiz * rm;

    R_direct  = sqrt(xh^2 + (zs - zr)^2);
    R_surf    = sqrt(xh^2 + (zs + zr)^2);
    R_bot     = sqrt(xh^2 + (2*H - zs - zr)^2);
    R_surfbot = sqrt(xh^2 + (2*H - zs + zr)^2);

    % Reflection coefficients
    Gamma_surf = -1;

    grazing_angle_bot = atan2((2*H - zs - zr), xh); % approx
    arg_sqrt = (c1/c2)^2 - cos(grazing_angle_bot)^2;
    if arg_sqrt >= 0
        Gamma_bot = ((rho2/rho1)*sin(grazing_angle_bot) - sqrt(arg_sqrt)) / ...
                     ((rho2/rho1)*sin(grazing_angle_bot) + sqrt(arg_sqrt));
    else
        Gamma_bot = ((rho2/rho1)*sin(grazing_angle_bot) - 1j*sqrt(-arg_sqrt)) / ...
                     ((rho2/rho1)*sin(grazing_angle_bot) + 1j*sqrt(-arg_sqrt));
    end

    paths = [ ...
        struct('R',R_direct,  'Gamma',1), ...
        struct('R',R_surf,    'Gamma',Gamma_surf), ...
        struct('R',R_bot,     'Gamma',Gamma_bot), ...
        struct('R',R_surfbot, 'Gamma',Gamma_surf*Gamma_bot)];

    x_m = zeros(Nt,1,'like',x_tx_col);

    eps0 = eps('double');
    for pidx = 1:numel(paths)
        Ri = max(paths(pidx).R,1.0);
        Gamma_i = paths(pidx).Gamma;

        % Delay and attenuation (monostatic: round-trip over this path)
        tau_i = 2*Ri / c;

        % Two-way propagation loss, consistent with tau_i and with the
        % project's canonical convention in transmission_loss.m
        % (TL_2way = 40*log10(R) + 2*alpha*(R/1000)).
        TL_dB = 40*log10(max(Ri,1)) + 2*alpha_dB_km*(Ri/1000);

        Ai = Gamma_i .* 10.^(-TL_dB/20);

        % Fractional delay using spline on real/imag (vectorized over time)
        tq = T - tau_i;
        % indices for spline mapping (treat t as uniformly sampled defensively)
        % If t is uniform, spline against sample points.
        xqRe = spline(T, real(x_tx_col), tq);
        xqIm = spline(T, imag(x_tx_col), tq);
        x_m = x_m + Ai .* (xqRe + 1j*xqIm);
    end

if any(~isfinite(x_m))
    error('multipath_model: Numerical instability detected.');
end

    x_n_multipath(:,m) = amp .* x_m;
end

end

