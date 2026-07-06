function TS_total = target_scattering(theta_asp_vec, Lt, a, a_cap, fc, c)
%TARGET_SCATTERING Aspect-dependent TS(theta_asp) with MATLAB sinc convention.
%
%   TS_total = target_scattering(theta_asp_vec, Lt, a, a_cap, fc, c)
%
% Inputs:
%   theta_asp_vec - aspect angle vector (radians), any shape.
%   Lt            - characteristic length (m), scalar >0.
%   a             - radius/characteristic radius (m), scalar >0.
%   a_cap         - end-cap radius/characteristic (m), scalar >0.
%   fc            - center frequency (Hz), scalar >0.
%   c             - sound speed (m/s), scalar >0.
%
% Output:
%   TS_total - target strength in dB, same shape as theta_asp_vec.

validateattributes(theta_asp_vec, {'numeric'},{'real','finite','nonempty','vector'}, mfilename,'theta_asp_vec',1);
validateattributes(Lt, {'numeric'},{'real','finite','scalar','positive'}, mfilename,'Lt',2);
validateattributes(a, {'numeric'},{'real','finite','scalar','positive'}, mfilename,'a',3);
validateattributes(a_cap, {'numeric'},{'real','finite','scalar','positive'}, mfilename,'a_cap',4);
validateattributes(fc, {'numeric'},{'real','finite','scalar','positive'}, mfilename,'fc',5);
validateattributes(c, {'numeric'},{'real','finite','scalar','positive'}, mfilename,'c',6);

th = theta_asp_vec(:); % column
k = 2*pi*fc / c;
lambda_c = c / fc;

% Broadside term
TS_broadside = 10*log10( a * (Lt.^2) / (2*lambda_c) );

% sinc convention: MATLAB sinc(x)=sin(pi x)/(pi x)
% Need sinc_arg = (k Lt sin(theta))/pi
sinc_arg = (k .* Lt .* sin(th)) ./ pi;

diffraction = 20*log10( abs(sinc(sinc_arg)) + eps('double') );
aspect_falloff = 10*log10( abs(cos(th)) + eps('double') );

TS_env = TS_broadside + diffraction + aspect_falloff;

TS_end = 10*log10( (a_cap.^2)/4 + eps('double') );

phase_diff = 2*k*Lt*cos(th);
coherent_sum = 10.^(TS_env/20) + 10.^(TS_end/20) .* exp(1j*phase_diff);

TS_total_col = 10*log10( abs(coherent_sum).^2 + eps('double') );

TS_total = reshape(TS_total_col, size(theta_asp_vec));

end

