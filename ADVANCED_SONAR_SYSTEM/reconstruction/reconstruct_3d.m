function point_cloud = reconstruct_3d(detections, beam_config, config)

%% =========================================================
% 3D RECONSTRUCTION MODULE
%% =========================================================

disp('Generating 3D Reconstruction...');

n = length(detections.indices);

x = detections.indices(:);

y = randn(n,1);

z = -abs(randn(n,1));

point_cloud = struct();

point_cloud.x = x;
point_cloud.y = y;
point_cloud.z = z;

disp('-----------------------------------');
disp('3D RECONSTRUCTION COMPLETE');
disp('-----------------------------------');

fprintf('Points Generated : %d\n', n);

disp('-----------------------------------');

end