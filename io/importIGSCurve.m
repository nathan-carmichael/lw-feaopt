function result = importIGSCurve(filename, options)
% IMPORTIGSCURVE  Import a 2-D curve from an IGES / IGS file.
%
%   Reads the following IGES entity types:
%     110 — Line
%     100 — Circular Arc
%     104 — Conic Arc
%     112 — Parametric Spline
%     126 — B-Spline Curve  (most common for SolidWorks splines)
%     106 — Copious Data    (point sequences)
%
%   Entities are chained end-to-end, then resampled to uniform arc-length
%   spacing for FEA meshing.
%
%   result = importIGSCurve('sketch.igs')
%   result = importIGSCurve('sketch.igs', opts)
%
%   Options (struct, all optional)
%   ------------------------------
%   .num_elements        Number of FEA elements              (250)
%   .scale               Coordinate scale factor             (1)
%   .samples_per_entity  Fine sampling per entity            (300)
%   .flip                Reverse curve direction             (false)
%   .plot                Show diagnostic figure              (true)
%   .verbose             Print summary to console            (true)
%
%   Output
%   ------
%   result  struct — same layout as importCSVCurve:
%     .nodes, .connectivity, .num_nodes, .num_elements,
%     .arc_length, .raw_entities, .raw_points,
%     .start_tangent, .end_tangent, .start_angle, .end_angle
%
%   See also: importCSVCurve, StructureGeometry

% =====================================================================
%  Options
% =====================================================================
if nargin < 2, options = struct(); end

opt = applyDefaults(options, struct( ...
    'num_elements',       250,  ...
    'scale',              1,    ...
    'samples_per_entity', 300,  ...
    'flip',               false,...
    'plot',               true, ...
    'verbose',            true));

% =====================================================================
%  Read & classify IGES sections
% =====================================================================
all_lines = splitlines(fileread(filename));

d_lines = {};
p_lines = {};
for i = 1:numel(all_lines)
    ln = all_lines{i};
    if numel(ln) < 73, continue; end
    switch ln(73)
        case 'D', d_lines{end+1} = ln; %#ok<AGROW>
        case 'P', p_lines{end+1} = ln; %#ok<AGROW>
    end
end

% =====================================================================
%  Parse Directory Entries (D-section, two lines per entity)
% =====================================================================
num_de = floor(numel(d_lines) / 2);
dir_entries = struct('entity_type',{}, 'param_start',{}, 'param_count',{});

for i = 1:num_de
    l1 = d_lines{2*i - 1};
    l2 = d_lines{2*i};
    de.entity_type = str2double(strtrim(l1(1:8)));
    de.param_start = str2double(strtrim(l1(9:16)));
    pc = str2double(strtrim(l2(25:32)));
    if isnan(pc), pc = 1; end
    de.param_count = pc;
    dir_entries(end+1) = de; %#ok<AGROW>
end

% =====================================================================
%  Collect Parameter Data (P-section)
% =====================================================================
p_map = containers.Map('KeyType','int32','ValueType','char');
for i = 1:numel(p_lines)
    ln  = p_lines{i};
    seq = int32(str2double(strtrim(ln(74:80))));
    dat = strtrim(ln(1:64));
    if isKey(p_map, seq)
        p_map(seq) = [p_map(seq), dat];
    else
        p_map(seq) = dat;
    end
end

% =====================================================================
%  Build entity list (supported types only)
% =====================================================================
SUPPORTED = [100, 104, 106, 110, 112, 126];

entities = {};
for i = 1:numel(dir_entries)
    de = dir_entries(i);
    if ~ismember(de.entity_type, SUPPORTED), continue; end

    pstr = '';
    for p = de.param_start:(de.param_start + de.param_count - 1)
        if isKey(p_map, int32(p))
            pstr = [pstr, p_map(int32(p))]; %#ok<AGROW>
        end
    end

    pstr  = strrep(pstr, ';', '');
    parts = strsplit(pstr, ',');
    vals  = zeros(numel(parts), 1);
    for j = 1:numel(parts)
        v = str2double(strtrim(parts{j}));
        if isnan(v), v = 0; end
        vals(j) = v;
    end

    entities{end+1} = struct('type', de.entity_type, 'params', vals); %#ok<AGROW>
end

if opt.verbose
    fprintf('\n=== IGES Import: %s ===\n', filename);
    fprintf('  Directory entries: %d\n', numel(dir_entries));
    fprintf('  Curve entities:   %d\n', numel(entities));
    for i = 1:numel(entities)
        fprintf('    %d: type %d (%s)\n', i, entities{i}.type, ...
            igesTypeName(entities{i}.type));
    end
end

if isempty(entities)
    error('FEAOPT:noEntities', ...
        'No supported curve entities found in %s.', filename);
end

% =====================================================================
%  Sample each entity into a polyline
% =====================================================================
polys = cell(numel(entities), 1);
for i = 1:numel(entities)
    polys{i} = sampleEntity(entities{i}, opt.samples_per_entity);
end

empty = cellfun(@isempty, polys);
polys(empty)    = [];
entities(empty) = [];

if isempty(polys)
    error('FEAOPT:noData', 'No valid curve data extracted.');
end

% =====================================================================
%  Chain entities into a single polyline
% =====================================================================
[chained, chain_order] = chainEntities(polys, opt.verbose);

if isempty(chained)
    error('FEAOPT:chainFailed', 'Could not chain entities into a connected curve.');
end

chained = chained * opt.scale;

if opt.flip
    chained = flipud(chained);
end

% =====================================================================
%  Resample to uniform arc-length spacing
% =====================================================================
[nodes, total_arc] = resampleByArcLength(chained, opt.num_elements);

% =====================================================================
%  Connectivity & tangents
% =====================================================================
num_el    = opt.num_elements;
num_nodes = size(nodes, 1);
connectivity = [(1:num_el)', (2:num_el+1)'];

[start_tangent, start_angle] = edgeTangent(nodes(1,:),     nodes(2,:));
[end_tangent,   end_angle]   = edgeTangent(nodes(end-1,:), nodes(end,:));

% =====================================================================
%  Output struct
% =====================================================================
result.nodes         = nodes;
result.connectivity  = connectivity;
result.num_nodes     = num_nodes;
result.num_elements  = num_el;
result.arc_length    = total_arc;
result.raw_entities  = entities;
result.raw_points    = chained;
result.start_tangent = start_tangent;
result.end_tangent   = end_tangent;
result.start_angle   = start_angle;
result.end_angle     = end_angle;

if opt.verbose
    fprintf('\n  Arc length:    %.4f (scaled units)\n', total_arc);
    fprintf('  Nodes: %d  |  Elements: %d\n', num_nodes, num_el);
    fprintf('  Start: [%.4f, %.4f]  tang: %.1f deg\n', ...
        nodes(1,:), start_angle*180/pi);
    fprintf('  End:   [%.4f, %.4f]  tang: %.1f deg\n', ...
        nodes(end,:), end_angle*180/pi);
end

% =====================================================================
%  Diagnostic plot
% =====================================================================
if opt.plot
    plotIGSImportResult(polys, opt.scale, chained, nodes, total_arc, ...
        num_el, start_tangent, end_tangent);
end

end


%% ========================= IGES ENTITY SAMPLING =========================

function pts = sampleEntity(ent, n)
    switch ent.type
        case 110, pts = sampleLine110(ent.params, n);
        case 100, pts = sampleArc100(ent.params, n);
        case 126, pts = sampleBSpline126(ent.params, n);
        case 112, pts = sampleParamSpline112(ent.params, n);
        case 106, pts = sampleCopious106(ent.params);
        case 104, pts = sampleConic104(ent.params, n);
        otherwise, pts = [];
    end
end

function pts = sampleLine110(p, n)
    % Entity 110 — Line segment
    % p: [110, x1, y1, z1, x2, y2, z2]
    i = 2;
    pts = [linspace(p(i), p(i+3), n)', linspace(p(i+1), p(i+4), n)'];
end

function pts = sampleArc100(p, n)
    % Entity 100 — Circular arc (CCW)
    % p: [100, ZT, cx, cy, xs, ys, xe, ye]
    i  = 2;
    cx = p(i+1); cy = p(i+2);
    xs = p(i+3); ys = p(i+4);
    xe = p(i+5); ye = p(i+6);
    r  = sqrt((xs-cx)^2 + (ys-cy)^2);
    a0 = atan2(ys-cy, xs-cx);
    a1 = atan2(ye-cy, xe-cx);
    if a1 <= a0, a1 = a1 + 2*pi; end
    ang = linspace(a0, a1, n)';
    pts = [cx + r*cos(ang), cy + r*sin(ang)];
end

function pts = sampleBSpline126(p, n)
    % Entity 126 — Rational B-Spline (NURBS)
    i = 2;
    K = round(p(i));        % num_ctrl_pts - 1
    M = round(p(i+1));      % degree
    i = i + 6;              % skip PROP1..PROP4

    N = K + 1;              % number of control points
    A = N + M + 1;          % number of knots

    knots   = p(i : i+A-1);        i = i + A;
    weights = p(i : i+N-1);        i = i + N;

    cp = zeros(N, 2);
    for k = 1:N
        cp(k,:) = [p(i), p(i+1)];
        i = i + 3;                  % skip z
    end

    if i+1 <= numel(p)
        v0 = p(i); v1 = p(i+1);
    else
        v0 = knots(M+1); v1 = knots(end-M);
    end

    t_vals = linspace(v0, v1 - 1e-10, n)';
    pts = zeros(n, 2);
    for k = 1:n
        pts(k,:) = evalNURBS(t_vals(k), cp, knots, weights, M);
    end
end

function pt = evalNURBS(t, cp, knots, weights, deg)
    n_cp = size(cp, 1);

    % Find knot span
    k = deg + 1;
    while k < n_cp && knots(k+1) <= t
        k = k + 1;
    end

    % De Boor in homogeneous coordinates [w*x, w*y, w]
    d = zeros(deg+1, 3);
    for j = 0:deg
        ci = k - deg + j;
        if ci >= 1 && ci <= n_cp
            w = weights(ci);
            d(j+1,:) = [cp(ci,1)*w, cp(ci,2)*w, w];
        end
    end

    for r = 1:deg
        for j = deg:-1:r
            ci = k - deg + j;
            denom = knots(ci+deg-r+1) - knots(ci);
            if abs(denom) < 1e-14
                alpha = 0;
            else
                alpha = (t - knots(ci)) / denom;
            end
            d(j+1,:) = (1-alpha)*d(j,:) + alpha*d(j+1,:);
        end
    end

    hw = d(deg+1,:);
    if abs(hw(3)) > 1e-14
        pt = [hw(1)/hw(3), hw(2)/hw(3)];
    else
        pt = [hw(1), hw(2)];
    end
end

function pts = sampleParamSpline112(p, n)
    % Entity 112 — Parametric spline (cubic polynomial segments)
    i = 2;
    N_seg = round(p(i+3));
    i = i + 4;

    breaks = p(i : i+N_seg);
    i = i + N_seg + 1;

    pts = [];
    for s = 1:N_seg
        AX = p(i); BX = p(i+1); CX = p(i+2); DX = p(i+3); i = i+4;
        AY = p(i); BY = p(i+1); CY = p(i+2); DY = p(i+3); i = i+4;
        i = i + 4;  % skip Z coefficients

        dt = breaks(s+1) - breaks(s);
        tl = linspace(0, dt, max(2, round(n/N_seg)))';

        xseg = AX + BX*tl + CX*tl.^2 + DX*tl.^3;
        yseg = AY + BY*tl + CY*tl.^2 + DY*tl.^3;

        if isempty(pts)
            pts = [xseg, yseg];
        else
            pts = [pts; xseg(2:end), yseg(2:end)]; %#ok<AGROW>
        end
    end
end

function pts = sampleCopious106(p)
    % Entity 106 — Copious data (point list)
    i    = 2;
    flag = round(p(i));
    np   = round(p(i+1));
    i    = i + 2;

    pts = zeros(np, 2);
    switch flag
        case 1  % x,y pairs
            for k = 1:np, pts(k,:) = [p(i) p(i+1)]; i = i+2; end
        case 2  % x,y,z triples
            for k = 1:np, pts(k,:) = [p(i) p(i+1)]; i = i+3; end
        case 3  % x,y,z + vector
            for k = 1:np, pts(k,:) = [p(i) p(i+1)]; i = i+6; end
        otherwise
            pts = [];
    end
end

function pts = sampleConic104(p, n)
    % Entity 104 — Conic arc (simplified: straight line between endpoints)
    i = 2 + 6 + 1;   % skip coefficients + ZT
    x1 = p(i); y1 = p(i+1);
    x2 = p(i+2); y2 = p(i+3);
    pts = [linspace(x1,x2,n)', linspace(y1,y2,n)'];
end

function name = igesTypeName(t)
    map = containers.Map( ...
        {100,  104,         106,            110,    112,                  126}, ...
        {'Circular Arc','Conic Arc','Copious Data','Line','Parametric Spline','B-Spline'});
    if isKey(map, t), name = map(t); else, name = 'Unknown'; end
end


%% ========================= ENTITY CHAINING =========================

function [chained, order] = chainEntities(polys, verbose)
    n_ent = numel(polys);
    if n_ent == 0, chained = []; order = []; return; end
    if n_ent == 1, chained = polys{1}; order = 1; return; end

    starts = cell2mat(cellfun(@(p) p(1,:),   polys, 'Uniform', false));
    ends   = cell2mat(cellfun(@(p) p(end,:), polys, 'Uniform', false));

    used    = false(n_ent, 1);
    order   = zeros(n_ent, 1);
    flipped = false(n_ent, 1);

    order(1) = 1;  used(1) = true;
    cur_end  = ends(1,:);

    for step = 2:n_ent
        best_d = inf;  best_j = -1;  best_flip = false;

        for j = 1:n_ent
            if used(j), continue; end
            ds = norm(starts(j,:) - cur_end);
            de = norm(ends(j,:)   - cur_end);
            if ds < best_d, best_d = ds; best_j = j; best_flip = false; end
            if de < best_d, best_d = de; best_j = j; best_flip = true;  end
        end

        if best_j < 0
            warning('FEAOPT:chainGap', ...
                'Could not chain all entities. Using %d of %d.', step-1, n_ent);
            order = order(1:step-1);
            break
        end

        order(step) = best_j;
        used(best_j) = true;
        flipped(best_j) = best_flip;

        if best_flip
            cur_end = starts(best_j,:);
        else
            cur_end = ends(best_j,:);
        end

        if verbose && best_d > 1e-4
            fprintf('  Warning: gap of %.6f between entity %d and %d\n', ...
                best_d, order(step-1), best_j);
        end
    end

    chained = [];
    for step = 1:numel(order)
        pts = polys{order(step)};
        if flipped(order(step)), pts = flipud(pts); end
        if isempty(chained)
            chained = pts;
        else
            chained = [chained; pts(2:end,:)]; %#ok<AGROW>
        end
    end

    if verbose
        fprintf('  Chain order: ');
        for step = 1:numel(order)
            if flipped(order(step))
                fprintf('%d(flip) ', order(step));
            else
                fprintf('%d ', order(step));
            end
        end
        fprintf('\n');
    end
end


%% ========================= SHARED HELPERS =========================

function [nodes, total_arc] = resampleByArcLength(pts, num_el)
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
    d     = p2 - p1;
    tang  = d / norm(d);
    angle = atan2(tang(2), tang(1));
end

function plotIGSImportResult(polys, scale, chained, nodes, arc_len, ...
        num_el, t_start, t_end)
    figure('Name', 'IGES Import', 'Position', [100 100 800 400]);
    arrow_len = arc_len * 0.08;

    subplot(1,2,1); hold on;
    colors = hsv(numel(polys));
    for i = 1:numel(polys)
        pts = polys{i} * scale;
        plot(pts(:,1), pts(:,2), '-', 'LineWidth', 2, 'Color', colors(i,:));
    end
    plot(chained(1,1), chained(1,2), 'go', 'MarkerSize', 10, 'LineWidth', 2);
    plot(chained(end,1), chained(end,2), 'rs', 'MarkerSize', 10, 'LineWidth', 2);
    axis equal; grid on; xlabel('X'); ylabel('Y');
    title(sprintf('Raw Entities (%d)', numel(polys)));

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