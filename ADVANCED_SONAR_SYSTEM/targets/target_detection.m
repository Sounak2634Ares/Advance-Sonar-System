function detections = target_detection(processed_data, config)
% TARGET_DETECTION Vectorized 2D CA-CFAR target detection engine.
%
% Inputs:
%   processed_data - [NumSamples x NumBeams] numeric real energy matrix
%   config         - System configuration structure
%
% Outputs:
%   detections     - [M x 3] Matrix layout matching dashboard requirements:
%                    [Range_meters, Bearing_radians, NormalizedScore]

    % =========================================================================
    % 1. CRITICAL INITIALIZATION (Guarantees output variable always exists)
    % =========================================================================
    % This line ensures that even if 0 targets are detected, a clean empty 
    % matrix is returned to main.m instead of an unassigned variable error.
    detections = zeros(0, 3); 

    if nargin < 1 || isempty(processed_data)
        return;
    end

    [NumSamples, NumBeams] = size(processed_data);

    % Unpack configurations or apply safe defaults
    if isfield(config, 'Fs'),  Fs = config.Fs;   else, Fs = 48e3;  end
    if isfield(config, 'c'),   c = config.c;     else, c = 1500;  end
    if isfield(config, 'pfa'), Pfa = config.pfa; else, Pfa = 1e-4; end

    % Set up sliding window cell constraints
    guard_cells_r = 2;  guard_cells_b = 1;
    train_cells_r = 5;  train_cells_b = 2;

    % Calculate threshold multiplier factor
    N_train = (2*train_cells_r + 2*guard_cells_r + 1) * (2*train_cells_b + 2*guard_cells_b + 1) ...
              - (2*guard_cells_r + 1) * (2*guard_cells_b + 1);
    alpha_cfar = N_train * (Pfa^(-1/N_train) - 1);

    % Allow manual threshold overrides from configuration layer if specified
    if isfield(config, 'cfar_threshold_dB')
        alpha_cfar = 10^(config.cfar_threshold_dB / 10);
    end

    % =========================================================================
    % 2. 2D CA-CFAR CONVOLUTION SEARCH
    % =========================================================================
    kernel_full = ones(2*(train_cells_r + guard_cells_r) + 1, 2*(train_cells_b + guard_cells_b) + 1);
    kernel_train = kernel_full;
    
    r_start = train_cells_r + 1;
    r_end   = r_start + 2*guard_cells_r;
    b_start = train_cells_b + 1;
    b_end   = b_start + 2*guard_cells_b;
    kernel_train(r_start:r_end, b_start:b_end) = 0;
    kernel_train = kernel_train / sum(kernel_train(:));

    % Estimate background spatial noise floor map via high-speed convolution
    noise_floor = conv2(processed_data, kernel_train, 'same');
    
    % Core logical mask generation
    detection_mask = processed_data > (noise_floor * alpha_cfar);

    % Mute window boundaries to block convolution edge artifacts
    edge_r = train_cells_r + guard_cells_r;
    edge_b = train_cells_b + guard_cells_b;
    detection_mask(1:edge_r, :) = 0;
    detection_mask(end-edge_r:end, :) = 0;
    detection_mask(:, 1:edge_b) = 0;
    detection_mask(:, end-edge_b:end) = 0;

    % Find coordinates of valid targets
    [row_indices, col_indices] = find(detection_mask);

    % =========================================================================
    % 3. POPULATE COMPATIBLE MATRIX OUTPUT
    % =========================================================================
    if ~isempty(row_indices)
        % Map indices to physical parameters (Range and Bearing Grid)
        theta_grid = linspace(-pi/3, pi/3, NumBeams);
        ranges_vec = (0:NumSamples-1).' * (c / (2 * Fs));

        % Extract raw signal intensity scores
        raw_scores = processed_data(detection_mask);
        
        % Normalize scores safely between [0.0, 1.0] for the UI slider limits
        max_s = max(raw_scores);
        min_s = min(raw_scores);
        if max_s > min_s
            norm_scores = (raw_scores - min_s) / (max_s - min_s);
        else
            norm_scores = ones(size(raw_scores));
        end

        % Construct matrix layout: Column 1 = Range, Column 2 = Bearing, Column 3 = Score
        detections = [ranges_vec(row_indices), theta_grid(col_indices).', norm_scores];
    end
end