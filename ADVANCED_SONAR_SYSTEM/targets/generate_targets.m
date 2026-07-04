function targets = generate_targets(config)

%% =========================================================
% TARGET GENERATION MODULE
%
% Purpose:
% Creates underwater targets for sonar detection.
%
% Current Targets:
% - Submarine
% - Rock
% - Mine
% - Shipwreck
%
% Future:
% - Fish schools
% - Moving targets
% - Doppler targets
% - AI classification labels
%% =========================================================

disp('Generating Underwater Targets...');

%% =========================================================
% LOAD PARAMETERS
%% =========================================================

max_range = config.max_range;
max_depth = config.max_depth;

%% =========================================================
% TARGET STRUCTURE
%% =========================================================

targets = struct();

%% =========================================================
% TARGET 1 : SUBMARINE
%% =========================================================

targets(1).id = 1;
targets(1).type = 'Submarine';

targets(1).x = 25;
targets(1).y = -10;
targets(1).z = -18;

targets(1).rcs = 1.0;
targets(1).size = 8;

%% =========================================================
% TARGET 2 : ROCK
%% =========================================================

targets(2).id = 2;
targets(2).type = 'Rock';

targets(2).x = 45;
targets(2).y = 15;
targets(2).z = -22;

targets(2).rcs = 0.6;
targets(2).size = 3;

%% =========================================================
% TARGET 3 : MINE
%% =========================================================

targets(3).id = 3;
targets(3).type = 'Mine';

targets(3).x = 65;
targets(3).y = -20;
targets(3).z = -15;

targets(3).rcs = 0.8;
targets(3).size = 2;

%% =========================================================
% TARGET 4 : SHIPWRECK
%% =========================================================

targets(4).id = 4;
targets(4).type = 'Shipwreck';

targets(4).x = 85;
targets(4).y = 10;
targets(4).z = -30;

targets(4).rcs = 1.5;
targets(4).size = 12;

%% =========================================================
% TARGET VALIDATION
%% =========================================================

for k = 1:length(targets)

    targets(k).x = min(max(targets(k).x,0),max_range);

    targets(k).z = max(targets(k).z,-max_depth);

end

%% =========================================================
% DISPLAY TARGET INFORMATION
%% =========================================================

disp('-----------------------------------');
disp('TARGETS GENERATED');
disp('-----------------------------------');

fprintf('Total Targets : %d\n', length(targets));

for k = 1:length(targets)

    fprintf('\n');

    fprintf('Target %d\n', targets(k).id);
    fprintf('Type : %s\n', targets(k).type);

    fprintf('Position : (%.1f , %.1f , %.1f)\n', ...
        targets(k).x, ...
        targets(k).y, ...
        targets(k).z);

end

disp('-----------------------------------');

end