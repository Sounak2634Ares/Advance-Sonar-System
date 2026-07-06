function gpu_status = initialize_gpu(config)

%% =========================================================
% GPU INITIALIZATION
%% =========================================================

disp('Checking GPU Availability...');

if config.use_gpu

    gpuDevice;

    gpu_status = true;

    disp('GPU Successfully Initialized.');

else

    gpu_status = false;

    disp('GPU Disabled.');

end

end