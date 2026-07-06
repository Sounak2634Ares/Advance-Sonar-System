function a_realistic = aperture_physics(theta_vec, N, d, c, fc, w_e, eta, gamma_c)
%APERTURE_PHYSICS Coupling-corrected, dilation-aware manifold vector library.
%
%   a_realistic = aperture_physics(theta_vec, N, d, c, fc, w_e, eta, gamma_c)
%
% Inputs:
%   theta_vec - angles in radians, vector.
%   N         - number of array elements (positive integer).
%   d         - element spacing (m), scalar >0.
%   c         - sound speed (m/s), scalar >0.
%   fc        - center frequency (Hz), scalar >0.
%   w_e       - effective aperture parameter (m), scalar >0.
%   eta       - coupling gain (scale), scalar.
%   gamma_c   - coupling exponential decay, scalar >=0.
%
% Output:
%   a_realistic - [N x numel(theta_vec)] manifold library.

validateattributes(theta_vec, {'numeric'},{'real','finite','vector','nonempty'}, mfilename,'theta_vec',1);
validateattributes(N, {'numeric'},{'real','finite','scalar','integer','positive'}, mfilename,'N',2);
validateattributes(d, {'numeric'},{'real','finite','scalar','positive'}, mfilename,'d',3);
validateattributes(c, {'numeric'},{'real','finite','scalar','positive'}, mfilename,'c',4);
validateattributes(fc, {'numeric'},{'real','finite','scalar','positive'}, mfilename,'fc',5);
validateattributes(w_e, {'numeric'},{'real','finite','scalar','positive'}, mfilename,'w_e',6);
validateattributes(eta, {'numeric'},{'real','finite','scalar'}, mfilename,'eta',7);
validateattributes(gamma_c, {'numeric'},{'real','finite','scalar','nonnegative'}, mfilename,'gamma_c',8);

theta_vec = theta_vec(:).';
K = numel(theta_vec);

lambda_c = c / fc;

% Mutual coupling Toeplitz matrix C (reuse across theta)
idx = 0:(N-1);
k_idx = idx(:).';
% Build first column for Toeplitz via k_idx
c_k = eta .* exp(-gamma_c.*idx(:)) .* exp(-1j*2*pi*fc.*idx(:).*d/c);
C = toeplitz(c_k, conj(c_k));

% Centered element index n_idx (column)
n_idx = (0:(N-1))' - (N-1)/2;

% Ideal steering
sinTh = sin(theta_vec);
% a_ideal: N x K
Aideal = exp(-1j*pi*n_idx * sinTh);

% Piston/directivity envelope E(theta) with MATLAB sinc argument convention
u = (w_e .* sinTh) / lambda_c; % 1 x K
E_theta = sinc(u);            % 1 x K

% Apply coupling: C * a_ideal
CA = C * Aideal;             % N x K

% Multiply by scalar envelope per column
a_realistic = CA .* E_theta;

end

