function app = launch_dashboard(environment)
%LAUNCH_DASHBOARD Sonar GUI dashboard with waterfall + 3D target visualization.
%
%   app = launch_dashboard(environment)
%
% The returned struct contains:
%   - app.ui: uifigure
%   - app.axesWaterfall: uiaxes
%   - app.axes3D: uiaxes
%   - app.thresholdSlider, app.gainSlider
%   - app.runOnceFcn: function handle to update visualization
%
% The runOnceFcn signature:
%   runOnceFcn(args)
% where args is a struct with fields:
%   - beam_energy   : [NumSamples x NumBeams]
%   - theta_vec     : [1 x NumBeams]
%   - range_vec     : [NumSamples x 1]
%   - detections    : detections array (format supported by reconstruct_3d)



    arguments
        environment struct
    end

    % ---------------------------
    % Build UI Dashboard Window Layout
    % ---------------------------
    app = struct();

    app.ui = uifigure('Name','Advanced Sonar Dashboard (Phase 3c)','Color',[0.03 0.04 0.05], ...
        'Position',[100 100 1200 700]);

    % Panels
    app.panelLeft = uipanel(app.ui, 'Title','Sonar Waterfall', ...
        'Position',[15 15 760 670], 'BackgroundColor',[0.03 0.04 0.05], 'ForegroundColor',[0.85 0.9 1]);

    app.panelRight = uipanel(app.ui, 'Title','Targets & Controls', ...
        'Position',[785 15 400 670], 'BackgroundColor',[0.03 0.04 0.05], 'ForegroundColor',[0.85 0.9 1]);

    % Axes
    app.axesWaterfall = uiaxes(app.panelLeft, 'Position',[10 45 740 615], 'Color',[0 0 0]);
    app.axesWaterfall.XColor = [0.6 0.8 1];
    app.axesWaterfall.YColor = [0.6 0.8 1];
    app.axesWaterfall.GridAlpha = 0.15;

    app.axes3D = uiaxes(app.panelRight, 'Position',[15 345 370 290], 'Color',[0 0 0]);
    app.axes3D.XColor = [0.6 0.8 1];
    app.axes3D.YColor = [0.6 0.8 1];
    app.axes3D.ZColor = [0.6 0.8 1];
    app.axes3D.GridAlpha = 0.15;
    view(app.axes3D, 35, 25);

    % Controls & Sliders
    app.lblThreshold = uilabel(app.panelRight, 'Position',[15 285 250 22], ...
        'Text','Detection threshold', 'FontColor',[0.85 0.9 1]);
    app.thresholdSlider = uislider(app.panelRight, ...
        'Position',[15 265 370 3], ...
        'Limits',[0 1], 'Value',0.15); % Default low sensitivity to catch initial synthetic targets

    app.lblGain = uilabel(app.panelRight, 'Position',[15 220 250 22], ...
        'Text','Gain multiplier', 'FontColor',[0.85 0.9 1]);
    app.gainSlider = uislider(app.panelRight, ...
        'Position',[15 200 370 3], ...
        'Limits',[0.25 4], 'Value',1);

    app.dynamicRangeLabel = uilabel(app.panelRight, 'Position',[15 145 250 22], ...
        'Text','Dynamic range (dB)', 'FontColor',[0.85 0.9 1]);
    app.dynamicRangeSlider = uislider(app.panelRight, ...
        'Position',[15 125 370 3], ...
        'Limits',[20 90], 'Value',55);

    app.lblLatency = uilabel(app.panelRight, 'Position',[15 80 200 22], ...
        'Text','Latency: -- ms', 'FontColor',[0.85 0.9 1]);
    app.lblTargetCount = uilabel(app.panelRight, 'Position',[15 55 200 22], ...
        'Text','Targets: --', 'FontColor',[0.85 0.9 1]);

    app.statusLabel = uilabel(app.panelRight, 'Position',[15 20 370 22], ...
        'Text','Status: idle', 'FontColor',[0.85 0.9 1]);

    % Render state mapping passed into render_visualization
    app.ui_state = struct();
    app.ui_state.polar_ax = app.axesWaterfall;
    app.ui_state.scatter_ax = app.axes3D;

    % Enable 3D interactions safely
    try
        app.axes3D.Interactions = ["rotate3d" "pan" "zoom"];
    catch
        % no-op for older MATLAB graphics versions
    end

    % ---------------------------
    % Real-time frame rendering engine
    % ---------------------------
    app.runOnceFcn = @runOnce;

    function runOnce(args)
        t0 = tic;

        if ~isfield(args,'beam_energy') || ~isfield(args,'theta_vec') || ~isfield(args,'detections')
            error('launch_dashboard.runOnceFcn: args must include beam_energy, theta_vec, detections');
        end

        beam_energy = args.beam_energy;
        theta_vec   = args.theta_vec;
        range_vec   = [];
        if isfield(args,'range_vec')
            range_vec = args.range_vec;
        end
        detections  = args.detections;

        % Robust dynamic slider threshold filtering
        thr = app.thresholdSlider.Value;
        detFiltered = detections;
        
        if ~isempty(detections) && size(detections, 2) >= 3
            scores = detections(:,3);
            detFiltered = detections(scores >= thr, :);
        elseif isempty(detections)
            detFiltered = double.empty(0, 3);
        end

        % Project coordinates safely into 3D Cartesian coordinates
        try
            detectionsXYZ = reconstruct_3d(detFiltered, environment, theta_vec, range_vec);
        catch
            try
                detectionsXYZ = reconstruct_3d(detFiltered, environment);
            catch
                detectionsXYZ = double.empty(0, 3);
            end
        end

        % Pull updated modifier states from dashboard parameters
        environment_render = environment;
        environment_render.gainDisplay = app.gainSlider.Value;
        environment_render.dynamicRange_dB = app.dynamicRangeSlider.Value;

        % Render matrices using your visualization framework
        app.ui_state = render_visualization(beam_energy, theta_vec, range_vec, detectionsXYZ, environment_render, app.ui_state);

        % Update monitoring stats on screen panel
        latencyMs = toc(t0) * 1000;
        app.lblLatency.Text = sprintf('Latency: %.1f ms', latencyMs);
        app.lblTargetCount.Text = sprintf('Targets: %d', size(detectionsXYZ,1));
        app.statusLabel.Text = 'Status: updated';

        drawnow limitrate;
        app.lastFrame = struct('latencyMs', latencyMs, 'targetCount', size(detectionsXYZ,1));
    end
end