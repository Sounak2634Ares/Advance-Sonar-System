function [v_n_matrix, R_vv] = ambient_noise_model(fc, wind_speed, shipping_factor, N, d, c, Ns)
%AMBIENT_NOISE_MODEL Wenz-curve composite ambient noise with spatial coherence.
%
%   [v_n_matrix, R_vv] = ambient_noise_model(fc, wind_speed, shipping_factor, N, d, c, Ns)
%
% Inputs:
%   fc             - center frequency (Hz)
%   wind_speed    - wind speed (m/s)
%   shipping_factor - factor in [0,1]
%   N              - number of array elements (positive int)
%   d              - element spacing (m)
%   c              - sound speed (m/s)
%   Ns             - number of time samples
%
% Outputs:
%   v_n_matrix - [Ns x N] correlated complex Gaussian noise
%   R_vv        - [N x N] noise covariance at fc

validateattributes(fc, {'numeric'},{'real','finite','scalar','positive'}, mfilename,'fc',1);
validateattributes(wind_speed, {'numeric'},{'real','finite','scalar','nonnegative'}, mfilename,'wind_speed',2);
validateattributes(shipping_factor, {'numeric'},{'real','finite','scalar'}, mfilename,'shipping_factor',3);
validateattributes(N, {'numeric'},{'real','finite','scalar','integer','positive'}, mfilename,'N',4);
validateattributes(d, {'numeric'},{'real','finite','scalar','positive'}, mfilename,'d',5);
validateattributes(c, {'numeric'},{'real','finite','scalar','positive'}, mfilename,'c',6);
validateattributes(Ns, {'numeric'},{'real','finite','scalar','integer','positive'}, mfilename,'Ns',7);

shipping_factor = min(max(shipping_factor,0),1);

f_khz = fc/1000;

NL_ship = 40 + 20*(shipping_factor-0.5) + 26*log10(f_khz) - 60*log10(f_khz+0.03);
NL_wind = 50 + 7.5*sqrt(wind_speed) + 20*log10(f_khz) - 40*log10(f_khz+0.4);
NL_thermal = -15 + 20*log10(f_khz);

% Turbulence term: negligible at typical fc; include defensively
if f_khz < 10/1000
    NL_turb = 17 - 30*log10(f_khz*1000); % approximate conversion; kept for completeness
else
    NL_turb = -Inf;
end

eps0 = eps('double');

p_ship = 10.^(NL_ship/10);
p_wind = 10.^(NL_wind/10);
p_therm = 10.^(NL_thermal/10);
if isfinite(NL_turb)
    p_turb = 10.^(NL_turb/10);
else
    p_turb = 0;
end

p_total = p_turb + p_ship + p_wind + p_therm + eps0;
NL_total_dB = 10*log10(p_total);

sigma_v2 = 10.^(NL_total_dB/10);

% Spatial coherence
idx = (0:N-1)';
diffIdx = abs(idx - idx');

u = (2*fc*d .* diffIdx) / c;
Gamma = sinc(u);

R_vv = sigma_v2 .* Gamma;

% Noise synthesis: correlated complex Gaussian
% Regularized Cholesky
L_chol = chol(R_vv + eps0*eye(N), 'lower');
noise_iid = (randn(Ns,N) + 1j*randn(Ns,N)) / sqrt(2);

v_n_matrix = noise_iid * L_chol';

end

