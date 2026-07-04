function render_visualization( ...
            processed_data, ...
            point_cloud, ...
            environment, ...
            config)

%% =========================================================
% VISUALIZATION MODULE
%% =========================================================

disp('Rendering Visualization...');

figure('Name','Processed Sonar Signal');

plot(abs(processed_data.compressed_signal));

xlabel('Sample');
ylabel('Amplitude');

grid on;

figure('Name','3D Point Cloud');

scatter3( ...
    point_cloud.x, ...
    point_cloud.y, ...
    point_cloud.z, ...
    20, ...
    abs(point_cloud.z), ...
    'filled');

xlabel('X');
ylabel('Y');
zlabel('Z');

grid on;

title('Sonar Point Cloud');

disp('Visualization Complete.');

end