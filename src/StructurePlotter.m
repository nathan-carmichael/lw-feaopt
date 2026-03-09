classdef StructurePlotter < handle
    % STRUCTUREPLOTTER  Visualization utilities for 2-D frame structures.
    %
    %   Elements are drawn as filled rectangles whose width matches the
    %   element thickness, giving a physically realistic view of the
    %   structure.
    %
    %   Core workflow:
    %     plotter = StructurePlotter(geometry);           % geometry only
    %     plotter = StructurePlotter(geometry, analysis); % + deformed shapes
    %
    %     plotter.plotStructure();                        % undeformed model
    %     plotter.plotDeformedShape();                    % after analysis
    %     plotter.plotStressDistribution('combined');     % stress contours
    %     plotter.plotLoads(load_case);                   % force arrows
    %
    %   The low-level helpers Plot_Rectangle and Plot_Beam are public so
    %   that custom scripts (e.g. the demo) can build bespoke figures.
    %
    %   See also: StructureGeometry, FEAAnalysis, LoadCase

    % =====================================================================
    %  References
    % =====================================================================
    properties (SetAccess = private)
        geometry        % StructureGeometry (required)
        analysis        % FEAAnalysis       (optional — needed for deformation / stress)
    end

    % =====================================================================
    %  Appearance settings
    % =====================================================================
    properties
        node_size       % marker size for nodes
        deform_scale    % deformation magnification factor
        colors          % struct of RGB triples (see initializeColors)
    end

    % =====================================================================
    %  Constructor
    % =====================================================================
    methods
        function obj = StructurePlotter(geometry, analysis)
            % STRUCTUREPLOTTER  Create a plotter for the given geometry.
            %
            %   plotter = StructurePlotter(geometry)
            %   plotter = StructurePlotter(geometry, analysis)

            if ~isa(geometry, 'StructureGeometry')
                error('FEAOPT:badInput', ...
                    'First argument must be a StructureGeometry object.');
            end

            obj.geometry = geometry;

            if nargin > 1 && isa(analysis, 'FEAAnalysis')
                obj.analysis = analysis;
            else
                obj.analysis = [];
            end

            obj.node_size    = 6;
            obj.deform_scale = 1;
            obj.initializeColors();
        end
    end

    % =====================================================================
    %  Low-level drawing primitives
    % =====================================================================
    methods
        function Plot_Rectangle(~, x, y, thickness, color)
            % PLOT_RECTANGLE  Draw a filled rectangle between two points.
            %
            %   plotter.Plot_Rectangle([x1 x2], [y1 y2], thickness, color)
            %
            %   The rectangle is centred on the line from (x1,y1) to
            %   (x2,y2) with the given thickness measured perpendicular
            %   to that line.

            theta = pi/2 - atan2(y(2)-y(1), x(2)-x(1));
            ht = thickness / 2;
            ct = cos(theta);
            st = sin(theta);

            X = [x(1)+ht*ct, x(1)-ht*ct, x(2)-ht*ct, x(2)+ht*ct];
            Y = [y(1)-ht*st, y(1)+ht*st, y(2)+ht*st, y(2)-ht*st];

            patch(X, Y, color, 'LineStyle', 'none');
        end

        function Plot_Beam(obj, x, y, thickness, color)
            % PLOT_BEAM  Draw one or more connected beam segments.
            %
            %   plotter.Plot_Beam(x, y, thickness, color)
            %
            %   x, y       — vectors of node coordinates (length N)
            %   thickness  — scalar or [N-1 x 1] per segment
            %   color      — [1 x 3] RGB  or  [N-1 x 3] per segment

            for k = 1:numel(y)-1
                if size(color,1) > 1
                    c = color(k,:);
                else
                    c = color;
                end

                if numel(thickness) > 1
                    t = thickness(k);
                else
                    t = thickness;
                end

                obj.Plot_Rectangle([x(k) x(k+1)], [y(k) y(k+1)], t, c);
            end
        end
    end

    % =====================================================================
    %  Support symbols
    % =====================================================================
    methods
        function plotSupports(obj)
            % PLOTSUPPORTS  Draw constraint symbols at supported nodes.
            %
            %   Fully fixed (all 3 DOFs) — filled triangle
            %   Pinned (x + y only)      — circle

            sz = 0.02 * max(range(obj.geometry.nodes));

            % Build per-node constraint map
            node_con = false(obj.geometry.num_nodes, 3);
            for i = 1:numel(obj.geometry.constraints)
                c = obj.geometry.constraints(i);
                node_con(c.node_id, c.dof) = true;
            end

            for n = 1:obj.geometry.num_nodes
                if ~any(node_con(n,:)), continue; end

                xy = obj.geometry.getNodeCoords(n);

                if all(node_con(n,:))
                    % Fixed — triangle
                    tx = xy(1) + [-1  1  0 -1] * sz;
                    ty = xy(2) + [-1 -1 0.5 -1] * sz;
                    fill(tx, ty, 'k');
                elseif node_con(n,1) && node_con(n,2)
                    % Pinned — circle
                    rectangle('Position', [xy(1)-sz/2, xy(2)-sz/2, sz, sz], ...
                              'Curvature', [1 1], ...
                              'FaceColor', 'none', ...
                              'EdgeColor', 'k', 'LineWidth', 2);
                end
            end
        end
    end

    % =====================================================================
    %  Structure (undeformed)
    % =====================================================================
    methods
        function plotStructure(obj, opts)
            % PLOTSTRUCTURE  Draw the undeformed model with realistic thickness.
            %
            %   plotter.plotStructure()
            %   plotter.plotStructure(opts)
            %
            %   Options (struct):
            %     .show_nodes          logical  (true)
            %     .show_node_numbers   logical  (false)
            %     .show_member_numbers logical  (false)
            %     .show_supports       logical  (true)
            %     .title               string   ('Structure Geometry')

            if nargin < 2, opts = struct(); end
            opts = applyDefaults(opts, struct( ...
                'show_nodes', true, ...
                'show_node_numbers', false, ...
                'show_member_numbers', false, ...
                'show_supports', true, ...
                'title', 'Structure Geometry'));

            figure('Name', 'Structure Plot', 'Position', [100 100 800 600]);
            hold on;

            geom = obj.geometry;

            % Draw elements
            for i = 1:geom.num_members
                [n1, n2] = geom.getMemberNodes(i);
                xc = [geom.nodes(n1,1), geom.nodes(n2,1)];
                yc = [geom.nodes(n1,2), geom.nodes(n2,2)];

                if geom.is_rigid(i)
                    c = obj.colors.rigid;
                else
                    c = obj.colors.member_original;
                end
                obj.Plot_Beam(xc, yc, geom.thickness(i), c);
            end

            % Element labels
            if opts.show_member_numbers
                for i = 1:geom.num_members
                    [n1, n2] = geom.getMemberNodes(i);
                    mx = mean(geom.nodes([n1 n2], 1));
                    my = mean(geom.nodes([n1 n2], 2));
                    text(mx, my, num2str(i), 'HorizontalAlignment', 'center', ...
                        'BackgroundColor', 'w', 'EdgeColor', 'k', 'FontSize', 8);
                end
            end

            % Nodes
            if opts.show_nodes
                plot(geom.nodes(:,1), geom.nodes(:,2), 'o', ...
                    'Color', obj.colors.node, 'MarkerSize', obj.node_size, ...
                    'MarkerFaceColor', obj.colors.node);
                if opts.show_node_numbers
                    for i = 1:geom.num_nodes
                        text(geom.nodes(i,1), geom.nodes(i,2), ...
                            ['  ' num2str(i)], 'FontSize', 9);
                    end
                end
            end

            if opts.show_supports, obj.plotSupports(); end

            axis equal; grid on;
            xlabel('X [m]'); ylabel('Y [m]');
            title(opts.title, 'FontSize', 12, 'FontWeight', 'bold');
            hold off;
        end
    end

    % =====================================================================
    %  Deformed shape
    % =====================================================================
    methods
        function plotDeformedShape(obj, scale, opts)
            % PLOTDEFORMEDSHAPE  Draw original (ghost) and deformed structure.
            %
            %   plotter.plotDeformedShape()            % auto-scaled
            %   plotter.plotDeformedShape(10)           % custom scale
            %   plotter.plotDeformedShape([], opts)     % with options
            %
            %   Options (struct):
            %     .show_original  logical  (true)
            %     .show_nodes     logical  (true)
            %     .title          string   ('Deformed Shape')

            obj.requireAnalysis();

            if nargin < 2 || isempty(scale)
                max_d = obj.analysis.current_state.getMaxDisplacement();
                span  = max(range(obj.geometry.nodes));
                scale = 0.1 * span / max(max_d, eps);
            end

            if nargin < 3, opts = struct(); end
            opts = applyDefaults(opts, struct( ...
                'show_original', true, ...
                'show_nodes', true, ...
                'title', 'Deformed Shape'));

            figure('Name', 'Deformed Shape', 'Position', [150 150 800 600]);
            hold on;

            geom = obj.geometry;
            def  = obj.analysis.getDeformedCoordinates(scale);

            % Ghost of original
            if opts.show_original
                for i = 1:geom.num_members
                    [n1, n2] = geom.getMemberNodes(i);
                    obj.Plot_Beam( ...
                        [geom.nodes(n1,1), geom.nodes(n2,1)], ...
                        [geom.nodes(n1,2), geom.nodes(n2,2)], ...
                        geom.thickness(i), obj.colors.member_ghost);
                end
            end

            % Deformed
            for i = 1:geom.num_members
                [n1, n2] = geom.getMemberNodes(i);
                if geom.is_rigid(i)
                    c = obj.colors.rigid;
                else
                    c = obj.colors.member_deformed;
                end
                obj.Plot_Beam( ...
                    [def(n1,1), def(n2,1)], ...
                    [def(n1,2), def(n2,2)], ...
                    geom.thickness(i), c);
            end

            if opts.show_nodes
                if opts.show_original
                    plot(geom.nodes(:,1), geom.nodes(:,2), 'o', ...
                        'Color', obj.colors.member_ghost, ...
                        'MarkerSize', obj.node_size * 0.7);
                end
                plot(def(:,1), def(:,2), 'o', ...
                    'Color', obj.colors.member_deformed, ...
                    'MarkerSize', obj.node_size, ...
                    'MarkerFaceColor', obj.colors.member_deformed);
            end

            obj.plotSupports();

            axis equal; grid on;
            xlabel('X [m]'); ylabel('Y [m]');
            title(sprintf('%s  (scale: %.1fx)', opts.title, scale), ...
                'FontSize', 12, 'FontWeight', 'bold');
            if opts.show_original
                legend('Original', 'Deformed', 'Location', 'best');
            end
            hold off;
        end
    end

    % =====================================================================
    %  Load arrows
    % =====================================================================
    methods
        function plotLoads(obj, load_case, scale)
            % PLOTLOADS  Draw force and moment arrows on the structure.
            %
            %   plotter.plotLoads(load_case)           % auto-scaled
            %   plotter.plotLoads(load_case, 0.005)    % custom scale

            if nargin < 3
                span     = max(range(obj.geometry.nodes));
                load_mag = load_case.getTotalLoad().force_magnitude;
                scale    = 0.1 * span / max(load_mag, 1);
            end

            hold on;

            for i = 1:numel(load_case.load_records)
                rec = load_case.load_records(i);
                xy  = obj.geometry.getNodeCoords(rec.node_id);

                fx = rec.fx * load_case.lambda;
                fy = rec.fy * load_case.lambda;
                mz = rec.mz * load_case.lambda;

                if abs(fx) > eps
                    quiver(xy(1), xy(2), fx*scale, 0, 0, ...
                        'Color', obj.colors.load, 'LineWidth', 2, 'MaxHeadSize', 0.5);
                end
                if abs(fy) > eps
                    quiver(xy(1), xy(2), 0, fy*scale, 0, ...
                        'Color', obj.colors.load, 'LineWidth', 2, 'MaxHeadSize', 0.5);
                end

                % Moment — small arc with arrowhead
                if abs(mz) > eps
                    arc_r = 0.5 * scale;
                    if mz > 0
                        th = linspace(0, pi/2, 20);
                    else
                        th = linspace(pi, 3*pi/2, 20);
                    end
                    plot(xy(1) + arc_r*cos(th), xy(2) + arc_r*sin(th), ...
                        'Color', obj.colors.load, 'LineWidth', 2);
                    plot(xy(1) + arc_r*cos(th(end)), xy(2) + arc_r*sin(th(end)), ...
                        '>', 'Color', obj.colors.load, 'MarkerSize', 8, ...
                        'MarkerFaceColor', obj.colors.load);
                end
            end

            hold off;
        end
    end

    % =====================================================================
    %  Stress distribution
    % =====================================================================
    methods
        function plotStressDistribution(obj, stress_type)
            % PLOTSTRESSDISTRIBUTION  Colour elements by stress magnitude.
            %
            %   plotter.plotStressDistribution('combined')
            %   plotter.plotStressDistribution('axial')
            %   plotter.plotStressDistribution('bending')
            %   plotter.plotStressDistribution('shear')
            %
            %   Rigid elements are drawn in gray.

            obj.requireAnalysis();

            if isempty(obj.analysis.member_stresses)
                obj.analysis.performStressAnalysis();
            end

            stresses = obj.analysis.member_stresses;

            switch lower(stress_type)
                case 'axial',    data = stresses.axial;     lbl = 'Axial Stress [MPa]';
                case 'bending',  data = stresses.bending;   lbl = 'Bending Stress [MPa]';
                case 'shear',    data = stresses.shear;     lbl = 'Shear Stress [MPa]';
                case 'combined', data = stresses.combined;  lbl = 'Combined Stress [MPa]';
                otherwise
                    error('FEAOPT:badOption', ...
                        'Unknown stress type ''%s''. Use axial|bending|shear|combined.', ...
                        stress_type);
            end

            figure('Name', 'Stress Distribution', 'Position', [200 200 900 600]);
            hold on;

            geom = obj.geometry;
            cmap = jet(256);
            colormap(cmap);

            max_s = max(abs(data(:))) / 1e6;   % MPa

            for i = 1:geom.num_members
                [n1, n2] = geom.getMemberNodes(i);
                xc = [geom.nodes(n1,1), geom.nodes(n2,1)];
                yc = [geom.nodes(n1,2), geom.nodes(n2,2)];

                if geom.is_rigid(i)
                    c = obj.colors.rigid;
                else
                    s_mpa = max(abs(data(i,:))) / 1e6;
                    cidx  = round(min(s_mpa / max(max_s, eps), 1) * 255) + 1;
                    cidx  = max(1, min(cidx, 256));
                    c     = cmap(cidx, :);
                end

                obj.Plot_Beam(xc, yc, geom.thickness(i), c);
            end

            cb = colorbar;
            cb.Label.String     = lbl;
            cb.Label.FontWeight = 'bold';
            caxis([0, max_s]);

            plot(geom.nodes(:,1), geom.nodes(:,2), 'ko', ...
                'MarkerSize', obj.node_size * 0.5, 'MarkerFaceColor', 'k');
            obj.plotSupports();

            axis equal; grid on;
            xlabel('X [m]'); ylabel('Y [m]');
            title(sprintf('Stress Distribution: %s', stress_type), ...
                'FontSize', 12, 'FontWeight', 'bold');
            hold off;
        end
    end

    % =====================================================================
    %  Private helpers
    % =====================================================================
    methods (Access = private)
        function initializeColors(obj)
            obj.colors.member_original = [0.2, 0.4, 0.8];   % blue
            obj.colors.member_deformed = [0.8, 0.2, 0.2];   % red
            obj.colors.member_ghost    = [0.85, 0.85, 0.85]; % light gray
            obj.colors.rigid           = [0.3, 0.3, 0.3];   % dark gray
            obj.colors.node            = [0.1, 0.1, 0.1];   % near black
            obj.colors.support         = [0, 0, 0];          % black
            obj.colors.load            = [0, 0.6, 0];        % green
        end

        function requireAnalysis(obj)
            if isempty(obj.analysis)
                error('FEAOPT:noAnalysis', ...
                    'This method requires an FEAAnalysis object (pass it to the constructor).');
            end
        end
    end
end

% =========================================================================
%  Module-level utility (outside classdef — accessible to all methods)
% =========================================================================
function s = applyDefaults(s, defaults)
    % Merge default fields into an options struct without overwriting
    % fields the caller already set.
    flds = fieldnames(defaults);
    for i = 1:numel(flds)
        if ~isfield(s, flds{i})
            s.(flds{i}) = defaults.(flds{i});
        end
    end
end