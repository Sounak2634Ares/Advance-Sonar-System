classdef SignalProcessingModuleTest < matlab.unittest.TestCase
    % SIGNALPROCESSINGMODULETEST Quality Assurance Unit Suite for Sonar Processing Pipeline
    
    methods (Test)
        
        function testOutputDimensions(testCase)
            % Test 1: Dimensions stability verification
            config.NumSamples = 2048;
            config.N = 8;
            config.Fs = 48000;
            config.alpha = 0.15;
            config.c = 1500;
            config.T = 0.01;
            
            sensor_data = randn(config.NumSamples, config.N) + 1j*randn(config.NumSamples, config.N);
            tx_signal = randn(128, 1) + 1j*randn(128, 1);
            
            [processed, ~] = signal_processing(sensor_data, tx_signal, config);
            
            testCase.verifySize(processed, [config.NumSamples, config.N]);
        end
        
        function testMatchedFilterPeak(testCase)
            % Test 2: Verify pulse compression location and mainlobe resolution
            config.NumSamples = 4096;
            config.N = 4;
            config.Fs = 96000;
            config.c = 1500;
            config.alpha = 0.11;
            config.T = 0.01;
            
            % Generate reference chirp
            t_pulse = (0:255).'/config.Fs;
            BW = 12000; 
            tx_signal = exp(1j * pi * (BW/config.T) * t_pulse.^2);
            
            % Insert target delay at exactly index 1000
            target_idx = 1000;
            sensor_data = zeros(config.NumSamples, config.N);
            sensor_data(target_idx:target_idx+255, :) = repmat(tx_signal, 1, config.N);
            
            [processed, ~] = signal_processing(sensor_data, tx_signal, config);
            [~, peak_index] = max(abs(processed(:, 1)));
            
            % Filter peaks at the termination of incoming coherent correlation block window
            expected_peak = target_idx + length(tx_signal) - 1;
            testCase.verifyEqual(peak_index, expected_peak, 'AbsTol', 2);
        end
        
        function testTVGZeroRangeBoundaries(testCase)
            % Test 3: Ensure Singularity Avoidance algorithms work smoothly at t=0
            config.NumSamples = 512;
            config.N = 2;
            config.Fs = 44100;
            config.c = 1500;
            config.alpha = 0.2;
            config.T = 0.005;
            
            sensor_data = zeros(config.NumSamples, config.N);
            tx_signal = randn(64, 1);
            
            [processed, ~] = signal_processing(sensor_data, tx_signal, config);
            
            testCase.verifyFalse(any(isnan(processed(:))));
            testCase.verifyFalse(any(isinf(processed(:))));
        end
        
    end
end