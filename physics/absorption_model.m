function alpha_vec = absorption_model(T, S, pH, D, c, f_vec)
%ABSORPTION_MODEL Francois-Garrison alpha(f) derivation packaged.
%
%   alpha_vec = absorption_model(T, S, pH, D, c, f_vec)
%
% Inputs:
%   T      - temperature (°C), scalar or array broadcastable with f_vec.
%   S      - salinity (ppt), scalar.
%   pH     - pH, scalar.
%   D      - depth (m), scalar.
%   c      - sound speed (m/s), scalar or compatible.
%   f_vec  - frequency vector in kHz (or scalar), size arbitrary.
%
% Output:
%   alpha_vec - absorption coefficient (dB/km), same size as f_vec.

validateattributes(T, {'numeric'},{'real','finite','nonempty','scalar'}, mfilename,'T',1);
validateattributes(S, {'numeric'},{'real','finite','scalar','nonnegative'}, mfilename,'S',2);
validateattributes(pH, {'numeric'},{'real','finite','scalar'}, mfilename,'pH',3);
validateattributes(D, {'numeric'},{'real','finite','scalar','nonnegative'}, mfilename,'D',4);
validateattributes(c, {'numeric'},{'real','finite','scalar','positive'}, mfilename,'c',5);
validateattributes(f_vec, {'numeric'},{'real','finite','nonempty','vector'}, mfilename,'f_vec',6);

f_khz = f_vec(:).';

% Guard frequencies defensively
if any(f_khz < 0)
    error('%s: f_vec must be nonnegative (kHz).', mfilename);
end

% Parameters as specified
f1 = 0.78*sqrt(S/35) * exp(T/26);
A1 = (8.86/c) * 10.^(0.78*pH - 5);

f2 = 42 * exp(T/17);
A2 = 21.44*(S/c) * (1 + 0.025*T);
P2 = 1 - 1.37e-4*D + 6.2e-9*D.^2;

if T <= 20
    A3 = 4.937e-4 - 2.59e-5*T + 9.11e-7*T.^2 - 1.50e-8*T.^3;
else
    A3 = 3.964e-4 - 1.146e-5*T + 1.45e-7*T.^2 - 6.5e-10*T.^3;
end
P3 = 1 - 3.83e-5*D + 4.9e-10*D.^2;

% Compute alpha (vectorized)
alpha = (A1.*f1.*(f_khz.^2)) ./ (f_khz.^2 + f1^2) ...
      + (A2.*P2.*f2.*(f_khz.^2)) ./ (f_khz.^2 + f2^2) ...
      + A3.*P3.*(f_khz.^2);

% Preserve original shape
alpha_vec = reshape(alpha, size(f_vec));

end

