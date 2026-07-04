function launch_dashboard( ...
            processed_data, ...
            detections, ...
            point_cloud, ...
            config)

%% =========================================================
% DASHBOARD MODULE
%% =========================================================

disp('Launching Dashboard...');

fprintf('\n');
fprintf('=============================\n');
fprintf(' SONAR DASHBOARD\n');
fprintf('=============================\n');

fprintf('Detections : %d\n', ...
        length(detections.indices));

fprintf('Point Cloud Points : %d\n', ...
        length(point_cloud.x));

fprintf('=============================\n');

end