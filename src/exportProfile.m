function profile = exportProfile(geometry, options)
% EXPORTPROFILE  Rasterise an optimised structure, extract boundary loops,
%   fit splines, and export SolidWorks-ready point files.
%
%   The approach is topology-agnostic:
%     1. Draw every active element as a filled rectangle onto a binary image
%     2. Extract boundary contours using image processing (bwboundaries)
%     3. Convert pixel contours back to world coordinates
%     4. Fit a smooth spline to each closed loop
%     5. Export each loop as a separate text file
%
%   This correctly handles beams, trusses, branching structures, and any
%   topology with holes or multiple disconnected regions.
%
%   profile = exportProfile(geometry)
%   profile = exportProfile(geometry, opts)
%
%   Options (struct, all optional)
%   ------------------------------
%   .num_spline_points      Points per loop for spline fit      (1000)
%   .smoothing              csaps smoothing factor (0–1)        (0.999999999)
%   .scale_to_mm            Coordinate multiplier for export    (1000)
%   .resolution             Pixels per metre for rasterisation  (50000)
%   .min_thickness_frac     Exclude elements thinner than this
%                           fraction of max thickness            (0)
%   .export_filename        Base name for files ([] = skip)     ([])
%   .plot                   Show diagnostic figure              (true)
%   .verbose                Print summary                       (true)
%
%   Output
%   ------
%   profile  struct with fields:
%     .loops          {L x 1} cell — each entry is a struct with:
%                       .x, .y         raw pixel-traced boundary (world coords, scaled)
%                       .xs, .ys       spline-smoothed boundary (scaled)
%                       .is_hole       true if this loop is an inner hole
%                       .arc_length    total perimeter (scaled)
%     .num_loops      number of boundary loops found
%     .image          the binary silhouette image
%     .x_range        [xmin xmax] of the image in world coords
%     .y_range        [ymin ymax]
%     .scale          scale_to_mm value used
%     .export_files   cell array of written filenames
%
%   See also: StructureGeometry, optimizeThickness

% =====================================================================
%  Options
% =====================================================================
if nargin < 2, options = struct(); end

opt = applyDefaults(options, struct( ...
    'num_spline_points',  1000,        ...
    'smoothing',          0.9999999, ...
    'scale_to_mm',        1000,        ...
    'resolution',         50000,       ...
    'min_thickness_frac', 0,           ...
    'export_filename',    [],          ...
    'plot',               true,        ...
    'verbose',            true));

scale = opt.scale_to_mm;

% =====================================================================
%  Filter thin elements
% =====================================================================
thickness = geometry.thickness;
is_rigid  = geometry.is_rigid;
active    = ~is_rigid;

if opt.min_thickness_frac > 0
    t_max    = max(thickness(active));
    t_thresh = t_max * opt.min_thickness_frac;
    active   = active & (thickness >= t_thresh);

    if opt.verbose
        n_removed = sum(~is_rigid) - sum(active);
        if n_removed > 0
            fprintf('  Filtered %d thin elements (t < %.3f mm)\n', ...
                n_removed, t_thresh * scale);
        end
    end
end

active_ids = find(active);

if isempty(active_ids)
    error('FEAOPT:noActiveElements', 'No elements remain after filtering.');
end

% =====================================================================
%  Build element rectangles in world coordinates
% =====================================================================
rects = cell(numel(active_ids), 1);

for k = 1:numel(active_ids)
    i = active_ids(k);
    [n1, n2] = geometry.getMemberNodes(i);
    x1 = geometry.nodes(n1, 1);  y1 = geometry.nodes(n1, 2);
    x2 = geometry.nodes(n2, 1);  y2 = geometry.nodes(n2, 2);

    theta = atan2(y2 - y1, x2 - x1);
    ht    = geometry.thickness(i) / 2;
    px    = -sin(theta) * ht;
    py    =  cos(theta) * ht;

    % Four corners: CCW starting from start-node outer side
    rects{k} = [x1+px, y1+py;
                x2+px, y2+py;
                x2-px, y2-py;
                x1-px, y1-py];
end

% =====================================================================
%  Determine image bounds with padding
% =====================================================================
all_corners = cell2mat(rects);
x_min = min(all_corners(:,1));  x_max = max(all_corners(:,1));
y_min = min(all_corners(:,2));  y_max = max(all_corners(:,2));

% Add padding of 5% on each side
pad_x = (x_max - x_min) * 0.05;
pad_y = (y_max - y_min) * 0.05;
pad   = max(pad_x, pad_y);

x_min = x_min - pad;  x_max = x_max + pad;
y_min = y_min - pad;  y_max = y_max + pad;

% =====================================================================
%  Rasterise: draw filled rectangles onto binary image
% =====================================================================
res    = opt.resolution;   % pixels per metre
img_w  = max(2, round((x_max - x_min) * res));
img_h  = max(2, round((y_max - y_min) * res));

% Cap image size to avoid memory issues
max_dim = 8000;
if img_w > max_dim || img_h > max_dim
    shrink = max_dim / max(img_w, img_h);
    img_w  = round(img_w * shrink);
    img_h  = round(img_h * shrink);
    res    = img_w / (x_max - x_min);
    if opt.verbose
        fprintf('  Image capped at %dx%d pixels (effective res: %.0f px/m)\n', ...
            img_w, img_h, res);
    end
end

bw = false(img_h, img_w);

for k = 1:numel(rects)
    corners = rects{k};

    % Convert world coords to pixel coords
    px_x = round((corners(:,1) - x_min) / (x_max - x_min) * (img_w - 1)) + 1;
    px_y = round((corners(:,2) - y_min) / (y_max - y_min) * (img_h - 1)) + 1;

    % Clamp to image bounds
    px_x = max(1, min(img_w, px_x));
    px_y = max(1, min(img_h, px_y));

    % Flip y for image coordinates (row 1 = top)
    px_r = img_h - px_y + 1;

    % Create filled polygon mask and OR into image
    elem_mask = poly2mask(double(px_x), double(px_r), img_h, img_w);
    bw = bw | elem_mask;
end

if opt.verbose
    fprintf('  Rasterised %d elements onto %d x %d image\n', ...
        numel(active_ids), img_w, img_h);
end

% =====================================================================
%  Morphological cleanup: fill tiny gaps between adjacent elements
% =====================================================================
se = strel('disk', max(1, round(res * max(thickness(active)) * 0.1)));
bw = imclose(bw, se);

% Fill interior holes smaller than 1% of total area (noise)
bw = bwareaopen(bw, round(0.001 * img_w * img_h));

% =====================================================================
%  Extract boundary contours
% =====================================================================
[boundaries, labels, num_regions, adj] = bwboundaries(bw, 8, 'noholes');

% Also get boundaries WITH holes to identify inner loops
[boundaries_full, ~] = bwboundaries(bw, 8);

% Classify: outer boundaries vs holes
% bwboundaries with default returns outer first, then holes
% We use the full version and check enclosure
loop_data = cell(numel(boundaries_full), 1);

for b = 1:numel(boundaries_full)
    bd = boundaries_full{b};

    % Convert pixel coords back to world coords
    % bd is [row, col] — row 1 is top of image
    world_x = (bd(:,2) - 1) / (img_w - 1) * (x_max - x_min) + x_min;
    world_y = (img_h - bd(:,1)) / (img_h - 1) * (y_max - y_min) + y_min;

    % Determine if this is a hole: check if it's enclosed by another boundary
    % Simple heuristic: if the centroid of this boundary is inside bw after
    % filling, it's an outer boundary; otherwise compare with filled version
    is_hole = false;
    if b > num_regions
        is_hole = true;
    end

    loop_data{b} = struct( ...
        'world_x', world_x, ...
        'world_y', world_y, ...
        'is_hole', is_hole, ...
        'num_points', numel(world_x));
end

if opt.verbose
    n_outer = sum(~cellfun(@(L) L.is_hole, loop_data));
    n_holes = sum( cellfun(@(L) L.is_hole, loop_data));
    fprintf('  Found %d boundary loops (%d outer, %d holes)\n', ...
        numel(loop_data), n_outer, n_holes);
end

% =====================================================================
%  Fit splines to each loop
% =====================================================================
profile.loops = cell(numel(loop_data), 1);

for b = 1:numel(loop_data)
    ld = loop_data{b};
    wx = ld.world_x;
    wy = ld.world_y;

    % Remove duplicate consecutive points
    dx = diff(wx);
    dy = diff(wy);
    ds = sqrt(dx.^2 + dy.^2);
    keep = [true; ds > 1e-12];
    wx = wx(keep);
    wy = wy(keep);

    % Ensure closed
    if norm([wx(end)-wx(1), wy(end)-wy(1)]) > 1e-10
        wx(end+1) = wx(1); %#ok<AGROW>
        wy(end+1) = wy(1); %#ok<AGROW>
    end

    % Arc-length parameterisation
    ds_seg = sqrt(diff(wx).^2 + diff(wy).^2);
    s = [0; cumsum(ds_seg)];
    total_perim = s(end);

    % Spline fit
    if numel(s) >= 4 && total_perim > 0
        s_fine = linspace(0, total_perim, opt.num_spline_points + 1)';
        s_fine = s_fine(1:end-1);  % remove duplicate closing point

        try
            xs = csaps(s, wx, opt.smoothing, s_fine);
            ys = csaps(s, wy, opt.smoothing, s_fine);
        catch
            xs = interp1(s, wx, s_fine, 'pchip');
            ys = interp1(s, wy, s_fine, 'pchip');
        end

        % Close
        xs(end+1) = xs(1);
        ys(end+1) = ys(1);
    else
        xs = wx;
        ys = wy;
    end

    % Store (apply scale)
    loop = struct();
    loop.x          = wx * scale;
    loop.y          = wy * scale;
    loop.xs         = xs * scale;
    loop.ys         = ys * scale;
    loop.is_hole    = ld.is_hole;
    loop.arc_length = total_perim * scale;
    loop.num_raw_points = numel(wx);

    profile.loops{b} = loop;
end

profile.num_loops  = numel(loop_data);
profile.image      = bw;
profile.x_range    = [x_min x_max];
profile.y_range    = [y_min y_max];
profile.scale      = scale;

% =====================================================================
%  Export
% =====================================================================
profile.export_files = {};
if ~isempty(opt.export_filename)
    profile.export_files = writeFiles(profile, opt);
end

% =====================================================================
%  Plot
% =====================================================================
if opt.plot
    plotProfile(profile, geometry, active_ids, opt);
end

end


%% ========================= FILE EXPORT =========================

function files = writeFiles(profile, opt)
% Write each boundary loop as a separate XYZ text file.

    base  = opt.export_filename;
    files = {};

    for b = 1:profile.num_loops
        loop = profile.loops{b};

        if loop.is_hole
            tag = sprintf('hole_%d', b);
        else
            tag = sprintf('loop_%d', b);
        end

        % Spline version
        f = sprintf('%s_%s.txt', base, tag);
        writeXYZ(f, loop.xs, loop.ys);
        files{end+1} = f; %#ok<AGROW>
    end

    % MATLAB data
    f = [base '.mat'];
    export_data = profile; %#ok<NASGU>
    save(f, 'export_data');
    files{end+1} = f;

    if opt.verbose
        fprintf('\n  Exported files:\n');
        for i = 1:numel(files)
            fprintf('    %s\n', files{i});
        end
        fprintf('\n  SolidWorks import:\n');
        fprintf('    Insert > Curve > Curve Through XYZ Points\n');
        fprintf('    Import each loop file as a separate spline.\n');
        fprintf('    Close each loop, then Extrude or Cut as needed.\n');
    end
end

function writeXYZ(filename, x, y)
    fid = fopen(filename, 'w');
    if fid == -1
        error('FEAOPT:fileWrite', 'Cannot open %s for writing.', filename);
    end
    for i = 1:numel(x)
        fprintf(fid, '%.6f %.6f 0.000000\n', x(i), y(i));
    end
    fclose(fid);
end


%% ========================= PLOTTING =========================

function plotProfile(profile, ~, ~, ~)

    figure('Name', 'Exported Profile', 'Position', [100 100 600 500]);
    hold on;

    blue = [0.3 0.55 0.85];

    for b = 1:profile.num_loops
        loop = profile.loops{b};

        if loop.is_hole
            fill(loop.xs, loop.ys, [1 1 1], 'EdgeColor', blue, ...
                'LineWidth', 1.5, 'LineStyle', '--');
        else
            fill(loop.xs, loop.ys, [0.75 0.85 0.95], 'EdgeColor', blue, ...
                'LineWidth', 2, 'FaceAlpha', 0.5);
        end
    end

    axis equal; grid on;
    xlabel('X [mm]'); ylabel('Y [mm]');
    title('Smoothed Spline Loops');

end


%% ========================= UTILITY =========================

function s = applyDefaults(s, defaults)
    flds = fieldnames(defaults);
    for i = 1:numel(flds)
        if ~isfield(s, flds{i})
            s.(flds{i}) = defaults.(flds{i});
        end
    end
end