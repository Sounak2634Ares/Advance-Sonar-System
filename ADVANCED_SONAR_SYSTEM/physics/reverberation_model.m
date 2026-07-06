function out = reverberation_model(range_m_vec, SL, Sv, Ss, Sb, T, c, alpha_dB_km, psi_solid_angle, beam_width, boundary_ranges)
%REVERBERATION_MODEL Production-grade reverberation model.
%
% Computes total reverberation from:
%   • Volume
%   • Surface
%   • Bottom
%
% Fully vectorized and numerically stable.

%% ---------------- Input Validation ----------------

validateattributes(range_m_vec,{'numeric'},...
    {'real','finite','nonempty','vector'},...
    mfilename,'range_m_vec',1);

validateattributes(SL,{'numeric'},...
    {'real','finite','scalar'},...
    mfilename,'SL',2);

validateattributes(Sv,{'numeric'},...
    {'real','finite','scalar'},...
    mfilename,'Sv',3);

validateattributes(Ss,{'numeric'},...
    {'real','finite','scalar'},...
    mfilename,'Ss',4);

validateattributes(Sb,{'numeric'},...
    {'real','finite','scalar'},...
    mfilename,'Sb',5);

validateattributes(T,{'numeric'},...
    {'real','finite','scalar','nonnegative'},...
    mfilename,'T',6);

validateattributes(c,{'numeric'},...
    {'real','finite','scalar','positive'},...
    mfilename,'c',7);

validateattributes(alpha_dB_km,{'numeric'},...
    {'real','finite','scalar','nonnegative'},...
    mfilename,'alpha_dB_km',8);

validateattributes(psi_solid_angle,{'numeric'},...
    {'real','finite','scalar','nonnegative'},...
    mfilename,'psi_solid_angle',9);

validateattributes(beam_width,{'numeric'},...
    {'real','finite','scalar','nonnegative'},...
    mfilename,'beam_width',10);

if nargin<11
    boundary_ranges=[];
end

if ~isempty(boundary_ranges)
    validateattributes(boundary_ranges,{'struct'},{});
end

%% ---------------- Numerical Guards ----------------

R = range_m_vec(:);

% 1 mm minimum range
R_safe = max(R,1e-3);

alpha = alpha_dB_km;

%% ---------------- Volume Reverberation ----------------

geom_vol = -20*log10(R_safe);

abs_loss = -2*alpha.*R_safe/1000;

vol_term = Sv + 10*log10(max(c*T*psi_solid_angle/2,eps));

RL_vol = SL + geom_vol + abs_loss + vol_term;

%% ---------------- Boundary Reverberation ----------------

dR = c*T/2;

A = max(R_safe*dR*beam_width,eps);

RL_surf = SL ...
          -40*log10(R_safe) ...
          + abs_loss ...
          + Ss ...
          + 10*log10(A);

RL_bot = SL ...
         -40*log10(R_safe) ...
         + abs_loss ...
         + Sb ...
         + 10*log10(A);

%% ---------------- Optional Gating ----------------

if ~isempty(boundary_ranges)

    if isfield(boundary_ranges,'surface_mask')

        ms = logical(boundary_ranges.surface_mask);

        if ~isequal(size(ms),size(R))
            error('surface_mask size mismatch.');
        end

        RL_surf(~ms) = -Inf;

    end

    if isfield(boundary_ranges,'bottom_mask')

        mb = logical(boundary_ranges.bottom_mask);

        if ~isequal(size(mb),size(R))
            error('bottom_mask size mismatch.');
        end

        RL_bot(~mb) = -Inf;

    end

end

%% ---------------- Combine Powers ----------------

Pvol = 10.^(RL_vol/10);

Psurf = 10.^(RL_surf/10);

Pbot = 10.^(RL_bot/10);

Ptotal = Pvol + Psurf + Pbot;

Ptotal(~isfinite(Ptotal)) = eps;

RL_total_dB = 10*log10(Ptotal);

RL_total_dB(~isfinite(RL_total_dB)) = -300;

RL_total_linear = 10.^(RL_total_dB/20);

RL_total_linear(~isfinite(RL_total_linear)) = 0;

%% ---------------- Output ----------------

out = struct;

out.RL_dB = reshape(RL_total_dB,size(range_m_vec));

out.RL_linear = reshape(RL_total_linear,size(range_m_vec));

end