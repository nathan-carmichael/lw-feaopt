function result = importCSVCurve(filename, options)
% IMPORTCSVCURVE  Import a 2-D curve from a CSV of (x, y) coordinates.
%
%   The file should have two columns.  A header row is auto-detected and
%   skipped.  Points are assumed to be ordered along the curve and are
%   resampled to uniform arc-length spacing for FEA meshing.
%
%   result = importCSVCurve('curve.csv')
%   result = importCSVCurve('curve.csv', opts)
%
%   Options (struct, all optional)
%   ------------------------------
%   .num_elements  Number of elements along the curve         (250)
%   .scale         Coordinate scale factor (e.g. 0.001 mm→m) (1)
%   .flip          Reverse curve direction                    (false)
%   .plot          Show diagnostic figure                     (true)
%   .verbose       Print summary to console                   (true)
%
%   Output
%   ------
%   result  struct with fields:
%     .nodes, .connectivity, .num_nodes, .num_elements,
%     .arc_length, .raw_points,
%     .start_tangent, .end_tangent, .start_angle, .end_angle
%
%   See also: importIGSCurve, StructureGeometry

% =====================================================================
%  Options
% =====================================================================
if nargin < 2, options = struct(); end

opt = applyDefaults(options, struct( ...
    'num_elements', 250, ...
    'scale',        1,   ...
    'flip',         false, ...
    'plot',         true, ...
    'verbose',      true));

% =====================================================================
%  Read CSV
% =====================================================================
raw = readmatrix(filename);

% Keep only first two columns
if size(raw, 2) > 2
    raw = raw(:, 1:2);
end

% Drop rows with NaN (headers / bad lines)
raw = raw(~any(isnan(raw), 2), :);

if size(raw, 1) < 2
    error('FEAOPT:tooFewPoints', ...
        'Need at least 2 points, found %d.', size(raw, 1));
end

raw_points = raw * opt.scale;

if opt.flip
    raw_points = flipud(raw_points);
end

if opt.verbose
    fprintf('\n=== CSV Import: %s ===\n', filename);
    fprintf('  Raw points: %d\n', size(raw_points, 1));
end

% =====================================================================
%  Remove duplicate consecutive points
% =====================================================================
seg_len = sqrt(sum(diff(raw_points).^2, 2));
keep    = [true; seg_len > 1e-12];
raw_points = raw_points(keep, :);

if opt.verbose && sum(~keep) > 0
    fprintf('  Removed %d duplicate points\n', sum(~keep));
end

% =====================================================================
%  Resample to uniform arc-length spacing
% =====================================================================
[nodes, total_arc] = resampleByArcLength(raw_points, opt.num_elements);

% =====================================================================
%  Connectivity & tangents
% =====================================================================
num_el    = opt.num_elements;
num_nodes = size(nodes, 1);
connectivity = [(1:num_el)', (2:num_el+1)'];

[start_tangent, start_angle] = edgeTangent(nodes(1,:),   nodes(2,:));
[end_tangent,   end_angle]   = edgeTangent(nodes(end-1,:), nodes(end,:));

% =====================================================================
%  Output struct
% =====================================================================
result.nodes         = nodes;
result.connectivity  = connectivity;
result.num_nodes     = num_nodes;
result.num_elements  = num_el;
result.arc_length    = total_arc;
result.raw_points    = raw_points;
result.start_tangent = start_tangent;
result.end_tangent   = end_tangent;
result.start_angle   = start_angle;
result.end_angle     = end_angle;

% =====================================================================
%  Console output
% =====================================================================
if opt.verbose
    fprintf('  Arc length:    %.4f (scaled units)\n', total_arc);
    fprintf('  Nodes: %d  |  Elements: %d\n', num_nodes, num_el);
    fprintf('  Start: [%.4f, %.4f]  tang: %.1f deg\n', ...
        nodes(1,:), start_angle * 180/pi);
    fprintf('  End:   [%.4f, %.4f]  tang: %.1f deg\n', ...
        nodes(end,:), end_angle * 180/pi);
end

% =====================================================================
%  Diagnostic plot
% =====================================================================
if opt.plot
    plotImportResult(raw_points, nodes, total_arc, num_el, ...
        start_tangent, end_tangent, 'CSV Import');
end

end

% =========================================================================
%  Shared helpers
% =========================================================================

function [nodes, total_arc] = resampleByArcLength(pts, num_el)
% Resample a polyline to num_el+1 nodes with uniform arc-length spacing.

    dxy      = diff(pts);
    ds       = sqrt(dxy(:,1).^2 + dxy(:,2).^2);
    cum_s    = [0; cumsum(ds)];
    total_arc = cum_s(end);

    s_target = linspace(0, total_arc, num_el + 1)';
    nodes    = zeros(num_el + 1, 2);
    nodes(1,:)   = pts(1,:);
    nodes(end,:) = pts(end,:);

    for i = 2:num_el
        lo = find(cum_s <= s_target(i), 1, 'last');
        hi = min(lo + 1, size(pts, 1));
        if lo == hi
            nodes(i,:) = pts(lo,:);
        else
            frac = (s_target(i) - cum_s(lo)) / (cum_s(hi) - cum_s(lo));
            nodes(i,:) = pts(lo,:) + frac * (pts(hi,:) - pts(lo,:));
        end
    end
end

function [tang, angle] = edgeTangent(p1, p2)
% Unit tangent and angle from p1 toward p2.
    d     = p2 - p1;
    tang  = d / norm(d);
    angle = atan2(tang(2), tang(1));
end

function plotImportResult(raw, nodes, arc_len, num_el, t_start, t_end, fig_name)
% Two-panel diagnostic figure: raw points and resampled mesh.

    figure('Name', fig_name, 'Position', [100 100 800 400]);
    arrow_len = arc_len * 0.08;

    subplot(1,2,1); hold on;
    plot(raw(:,1), raw(:,2), 'k.-', 'MarkerSize', 6);
    plot(raw(1,1), raw(1,2), 'go', 'MarkerSize', 10, 'LineWidth', 2);
    plot(raw(end,1), raw(end,2), 'rs', 'MarkerSize', 10, 'LineWidth', 2);
    axis equal; grid on; xlabel('X'); ylabel('Y');
    title(sprintf('Raw Points (%d)', size(raw,1)));

    subplot(1,2,2); hold on;
    plot(nodes(:,1), nodes(:,2), 'b.-', 'MarkerSize', 4);
    plot(nodes(1,1), nodes(1,2), 'go', 'MarkerSize', 10, 'LineWidth', 2);
    plot(nodes(end,1), nodes(end,2), 'rs', 'MarkerSize', 10, 'LineWidth', 2);
    quiver(nodes(1,1), nodes(1,2), t_start(1)*arrow_len, t_start(2)*arrow_len, 0, ...
        'g', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    quiver(nodes(end,1), nodes(end,2), t_end(1)*arrow_len, t_end(2)*arrow_len, 0, ...
        'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    axis equal; grid on; xlabel('X'); ylabel('Y');
    title(sprintf('Resampled (%d elements)', num_el));
end

function s = applyDefaults(s, defaults)
    flds = fieldnames(defaults);
    for i = 1:numel(flds)
        if ~isfield(s, flds{i})
            s.(flds{i}) = defaults.(flds{i});
        end
    end
end