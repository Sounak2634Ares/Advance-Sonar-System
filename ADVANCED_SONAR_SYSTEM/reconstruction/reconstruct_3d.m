function [pointsXYZ, info] = reconstruct_3d(detections, environment, theta_vec, range_vec)
% RECONSTRUCT_3D Projects polar sonar telemetry into 3D Cartesian space coordinates.
%
% File        : reconstruct_3d.m
% Purpose     : Phase 3c Precision Geometric Scene Mapping Engine
% Author      : Sounak Jana / Gemini Head Coordinator
%
% Inputs:
%   detections   - [M x 3] matrix containing [Range_m, Bearing_rad, Score]
%                  or [M x 4] matrix containing [RangeIdx, BeamIdx, Range_m, Bearing_rad]
%   environment  - Configuration struct detailing environmental bounds
%   theta_vec    - [1 x NumBeams] Vector of steering grid bearings (optional)
%   range_vec    - [NumSamples x 1] Vector of time-to-range parameters (optional)
%
% Outputs:
%   pointsXYZ    - [M x 3] Comprehensive Cartesian cloud matrix: [X_m, Y_m, Z_m]
%   info         - Execution statistics metadata wrapper

    if nargin < 2 || isempty(environment), environment = struct(); end

    % 1. Handle Zero Detections Edge Case Gracefully
    if isempty(detections)
        pointsXYZ = zeros(0, 3);
        info.count = 0;
        info.timestamp = datetime('now');
        return;
    end

    % 2. Parse Telemetry Formats Defensively
    if ismatrix(detections) && size(detections, 2) == 4
        % Precision production array format: [RangeIdx, BeamIdx, Range_m, Bearing_rad]
        range_m = detections(:, 3);
        bearing_rad = detections(:, 4);
    elseif ismatrix(detections) && size(detections, 2) >= 2
        % Standard format: [Range_m, Bearing_rad, Score]
        range_m = detections(:, 1);
        bearing_rad = detections(:, 2);
    else
        error('reconstruct_3d:InvalidFormat', ...
            'Input data structure layout error. Expected an array of explicit telemetry coordinates.');
    end

    % Force variables into uniform column vector alignments
    range_m = double(range_m(:));
    bearing_rad = double(bearing_rad(:));

    % =========================================================================
    % 3. CRITICAL AXIS FIX: Aligns with render_visualization polar projection
    % =========================================================================
    X = range_m .* cos(bearing_rad);
    Y = range_m .* sin(bearing_rad);

    % 4. Depth Profile Synthesis (Z-Axis Matrix Resolution)
    if isfield(environment, 'seabedZFcn') && isa(environment.seabedZFcn, 'function_handle')
        Z = environment.seabedZFcn(X, Y);
    elseif isfield(environment, 'seabedZ') && isa(environment.seabedZ, 'function_handle')
        Z = environment.seabedZ(X, Y);
    else
        % Fallback onto steady reference standard planes
        if isfield(environment, 'projectionZ'),      Z0 = environment.projectionZ;
        elseif isfield(environment, 'Z0'),            Z0 = environment.Z0;
        else,                                         Z0 = -15; % Default baseline depth floor
        end
        Z = Z0 + zeros(size(X), 'like', X);
    end

    % 5. Build Unified Consolidated Output Coordinates Matrix Space
    pointsXYZ = [X, Y, Z];

    % Pack tracking information wrapper
    info.count = size(pointsXYZ, 1);
    info.timestamp = datetime('now');
end