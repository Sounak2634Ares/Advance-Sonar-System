function processed_data = signal_processing(rx_data, tx_signal, config)

%% =========================================================
% SIGNAL PROCESSING MODULE
%
% Purpose:
% Performs matched filtering (pulse compression)
% on received sonar echoes.
%
%% =========================================================

disp('Processing Sonar Data...');

%% =========================================================
% LOAD SIGNALS
%% =========================================================

rx_signal = rx_data.signal;

%% =========================================================
% MATCHED FILTER
%% =========================================================

matched_filter = conj(fliplr(tx_signal));

%% =========================================================
% PULSE COMPRESSION
%% =========================================================

compressed_signal = conv( ...
                        rx_signal, ...
                        matched_filter, ...
                        'same');

%% =========================================================
% NORMALIZATION
%% =========================================================

compressed_signal = ...
    compressed_signal ./ ...
    max(abs(compressed_signal));

%% =========================================================
% OUTPUT STRUCTURE
%% =========================================================

processed_data = struct();

processed_data.raw_signal = rx_signal;

processed_data.compressed_signal = compressed_signal;

%% =========================================================
% DISPLAY INFO
%% =========================================================

disp('-----------------------------------');
disp('SIGNAL PROCESSING COMPLETE');
disp('-----------------------------------');

fprintf('Input Samples  : %d\n', ...
        length(rx_signal));

fprintf('Output Samples : %d\n', ...
        length(compressed_signal));

disp('-----------------------------------');

end