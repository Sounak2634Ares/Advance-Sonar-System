function s_rx = doppler_physics(s_tx_complex, Fs, R_dot, c, fc, tau0)
%DOPPLER_PHYSICS Apply Doppler time scaling to a complex IQ signal.
%
%   s_rx = doppler_physics(s_tx_complex, Fs, R_dot, c, fc, tau0)
%
% Inputs:
%   s_tx_complex : Complex transmit signal (row or column vector)
%   Fs           : Sampling frequency (Hz)
%   R_dot        : Relative radial velocity (m/s)
%   c            : Speed of sound (m/s)
%   fc           : Carrier frequency (Hz)
%   tau0         : Initial propagation delay (s)
%
% Output:
%   s_rx         : Doppler shifted signal (same size/orientation as input)

%% ---------------- Input Validation ----------------

validateattributes(s_tx_complex, {'numeric'}, ...
    {'vector','nonempty','finite'}, mfilename,'s_tx_complex',1);

validateattributes(Fs, {'numeric'}, ...
    {'real','finite','scalar','positive'}, mfilename,'Fs',2);

validateattributes(R_dot, {'numeric'}, ...
    {'real','finite','scalar'}, mfilename,'R_dot',3);

validateattributes(c, {'numeric'}, ...
    {'real','finite','scalar','positive'}, mfilename,'c',4);

validateattributes(fc, {'numeric'}, ...
    {'real','finite','scalar','positive'}, mfilename,'fc',5);

validateattributes(tau0, {'numeric'}, ...
    {'real','finite','scalar'}, mfilename,'tau0',6);

%% Preserve original orientation

isRow = isrow(s_tx_complex);

x = s_tx_complex(:);
N = numel(x);

%% Time axis (column vectors)

p_orig = (0:N-1).';
t_orig = p_orig / Fs;

%% Doppler scaling

beta = (c - R_dot) / (c + R_dot);

if ~isfinite(beta) || beta <= 0

    s_rx = zeros(size(s_tx_complex),'like',s_tx_complex);
    return

end

%% Zero Doppler shortcut

if abs(R_dot) < 1e-12 && abs(tau0) < eps

    s_rx = s_tx_complex;
    return

end

%% Doppler frequency

fD = -2 * R_dot * fc / c;

%% Query locations

t_query = beta .* (t_orig - tau0);

p_query = t_query * Fs;

%% Allocate output

s_rx_col = zeros(N,1,'like',x);

valid = (p_query >= 0) & (p_query <= (N-1));

if any(valid)

    idx = p_query(valid);

    re = interp1(p_orig, real(x), idx, 'spline');

    im = interp1(p_orig, imag(x), idx, 'spline');

    s_rx_col(valid) = complex(re,im);

end

%% Carrier phase correction

carrier = exp(1j * 2*pi * fD * t_orig);

s_rx_col = s_rx_col .* carrier;

%% Restore orientation

if isRow
    s_rx = s_rx_col.';
else
    s_rx = s_rx_col;
end

end