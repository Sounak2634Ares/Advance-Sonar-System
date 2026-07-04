function detections = target_detection(processed_data, config)
% TARGET_DETECTION Identify and isolate target echoes from processed sonar data.
%
% File        : target_detection.m
% Updated with defensive type checking for input data matrices

% =========================================================================
% STRUCT INTEGRITY GATEWAY
% =========================================================================
if isstruct(processed_data)
    % Fallback if an upstream module wrapped the signal inside an object
    if isfield(processed_data, 'compressed_signal')
        signal_matrix = processed_data.compressed_signal;
    elseif isfield(processed_data, 'signal')
        signal_matrix = processed_data.signal;
    else
        error('TargetDetection:InvalidStruct', 'Input struct missing recognizable signal matrix fields.');
    end
else
    % processed_data is already the raw [NumSamples x N] numeric matrix array
    signal_matrix = processed_data;
end

% Compute the analytic signal envelope magnitude safely from the extracted matrix
signal_envelope = abs(signal_matrix);

% =========================================================================
% REST OF YOUR DETECTION LOGIC (CFAR, Thresholding, etc.)
% =========================================================================
% Update the rest of your function variables to use 'signal_envelope' 
% or 'signal_matrix' instead of 'processed_data.compressed_signal'.