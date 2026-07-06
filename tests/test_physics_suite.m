%% ========================================================================
% test_physics_suite.m
% Advanced Multi-Beam Sonar System
% Comprehensive Physics Layer Verification
% MATLAB R2024b
%% ========================================================================

clear;
clc;

% Ensure sibling module folders are on path, regardless of current
% working directory or how this script was invoked (run vs F5 vs called
% from elsewhere).
this_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_dir);
addpath(fullfile(repo_root, 'physics'));

fprintf('\n=============================================\n');
fprintf(' ADVANCED SONAR PHYSICS TEST SUITE\n');
fprintf('=============================================\n\n');

global pass fail
pass = 0;
fail = 0;

run(@test_water_properties);
run(@test_sound_speed_profile);
run(@test_absorption_model);
run(@test_transmission_loss);
run(@test_target_scattering);
run(@test_doppler);
run(@test_ambient_noise);
run(@test_aperture);
run(@test_reverberation);
run(@test_multipath);

fprintf('\n=============================================\n');
fprintf('TOTAL PASSED : %d\n',pass);
fprintf('TOTAL FAILED : %d\n',fail);
fprintf('=============================================\n');

%% ------------------------------------------------------------------------
function run(f)

global pass fail

try
    f();
    fprintf('[PASS] %s\n',func2str(f));
    pass=pass+1;
catch ME
    fprintf('[FAIL] %s\n',func2str(f));
    fprintf('       %s\n',ME.message);
    fail=fail+1;
end

end

%% ------------------------------------------------------------------------
function test_water_properties()

z=(0:500:5000)';
T=20-0.003*z;
S=35*ones(size(z));

P=water_properties(T,S,z,45);

assert(all(isfinite(P.c)));
assert(all(isfinite(P.rho)));
assert(length(P.c)==length(z));

end

%% ------------------------------------------------------------------------
function test_sound_speed_profile()

z=(0:500:5000)';
T=20-0.003*z;
S=35*ones(size(z));

P=water_properties(T,S,z,45);

sp=sound_speed_profile(P,z);

assert(all(isfinite(sp.c_local)));
assert(sp.c_eff>1400);

end

%% ------------------------------------------------------------------------
function test_absorption_model()

alpha=absorption_model(15,35,8,5000,1500,[1 10 50]);

assert(all(alpha>=0));
assert(all(isfinite(alpha)));

end

%% ------------------------------------------------------------------------
function test_transmission_loss()

TL=transmission_loss([0 1 100 1000],0.04,1500);

assert(all(isfinite(TL.TL_dB)));
assert(TL.TL_dB(end)>TL.TL_dB(2));

end

%% ------------------------------------------------------------------------
function test_target_scattering()

th=deg2rad([-90 0 90]);

TS=target_scattering(th,5,0.5,0.5,30000,1500);

assert(all(isfinite(TS)));

end

%% ------------------------------------------------------------------------
function test_doppler()

Fs=96000;

t=(0:1023)'/Fs;

x=exp(1j*2*pi*5000*t);

y=doppler_physics(x,Fs,0,1500,30000,0);

assert(isfloat(y));
assert(~isreal(y));
assert(all(isfinite(real(y))));
assert(all(isfinite(imag(y))));
assert(norm(x-y)<1e-6);

end

%% ------------------------------------------------------------------------
function test_ambient_noise()

[X,R]=ambient_noise_model(30000,5,0.5,16,0.05,1500,1024);

assert(all(size(X)==[1024 16]));
assert(all(size(R)==[16 16]));
assert(all(isfinite(X(:))));

end

%% ------------------------------------------------------------------------
function test_aperture()

A=aperture_physics(linspace(-pi/2,pi/2,181),16,0.05,1500,30000,0.8,0.1,0.5);

assert(size(A,1)==16);
assert(size(A,2)==181);

end

%% ------------------------------------------------------------------------
function test_reverberation()

R=(0:100:5000)';

out=reverberation_model(R,220,-70,-30,-25,0.01,1500,0.04,0.5,1,[]);

assert(all(isfinite(out.RL_dB)));
assert(all(isfinite(out.RL_linear)));

end

%% ------------------------------------------------------------------------
function test_multipath()

t=(0:1/48000:0.01)';
x=exp(1j*2*pi*3000*t);

targets.amp=1;

Y=multipath_model( ...
    t,...
    x,...
    100,...
    20,...
    20,...
    100,...
    1000,...
    1500,...
    1800,...
    1700,...
    targets,...
    30000,...
    0.04,...
    1500);

assert(size(Y,1)==length(t));
assert(all(isfinite(real(Y(:)))));
assert(all(isfinite(imag(Y(:)))));

end