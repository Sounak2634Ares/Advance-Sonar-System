function detections = target_detection(processed_data, config)

%% =========================================================
% TARGET DETECTION MODULE
%% =========================================================

disp('Running Target Detection...');

signal = abs(processed_data.compressed_signal);

threshold = 0.3 * max(signal);

idx = find(signal > threshold);

detections = struct();

detections.indices = idx;
detections.amplitudes = signal(idx);

disp('-----------------------------------');
disp('TARGET DETECTION COMPLETE');
disp('-----------------------------------');

fprintf('Detections Found : %d\n', length(idx));

disp('-----------------------------------');

end