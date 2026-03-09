%% demo_truss_optimization.m
%  ShepherdLab LW FEA-OPT — Truss subdivision & thickness optimisation
%
%  Defines a Warren truss, subdivides every member into many smaller
%  frame elements, then runs thickness optimisation.  The result is an
%  organic, topology-optimised shape where material concentrates along
%  load paths and thins away elsewhere.
%
%  See also: subdivideTruss, optimizeThickness, exportProfile,
%            StructureGeometry, FEAAnalysis

clear; clc; close all;
set(0, 'DefaultAxesFontName', 'Calibri', 'DefaultTextFontName', 'Calibri');

%% ==================== TRUSS DEFINITION ====================

% Geometry
span   = 0.300;          % total span [m]  (300 mm)
height = 0.100;          % truss depth [m] (100 mm)
n_panels = 6;            % number of bottom-chord panels

% Material (steel)
E         = 200e9;       % Young's modulus  [Pa]
ultimate  = 400e6;       % UTS  [Pa]
rho       = 7850;        % density  [kg/m^3]

% Cross-section
beam_width        = 0.010;    % out-of-plane width  [m]  (10 mm)
initial_thickness = 0.015;    % starting thickness  [m]  (3 mm)

% Loading
Fy_load = 50000;         % vertical point load at centre bottom  [N]

% Subdivision
num_sub = 15;            % sub-elements per original member

%% ==================== BUILD TRUSS NODES & CONNECTIVITY ====================

% Bottom chord nodes  (n_panels + 1 nodes along y = 0)
dx_panel = span / n_panels;
n_bot    = n_panels + 1;
bot_x    = (0:n_panels)' * dx_panel;
bot_y    = zeros(n_bot, 1);

% Top chord nodes  (n_panels nodes, offset by half-panel, at y = height)
%   Warren truss: top nodes sit at panel midpoints
n_top  = n_panels;
top_x  = ((0:n_top-1)' + 0.5) * dx_panel;
top_y  = height * ones(n_top, 1);

% Combine:  bottom nodes first (1..n_bot), then top (n_bot+1..n_bot+n_top)
truss_nodes = [bot_x, bot_y;
               top_x, top_y];

n_total = size(truss_nodes, 1);
top_offset = n_bot;   % first top node index

% Build connectivity
conn = [];

% Bottom chord
for i = 1:n_bot-1
    conn(end+1, :) = [i, i+1]; %#ok<AGROW>
end

% Top chord
for i = 1:n_top-1
    conn(end+1, :) = [top_offset+i, top_offset+i+1]; %#ok<AGROW>
end

% Diagonals — Warren pattern: V shapes between bottom and top nodes
for i = 1:n_top
    % Rising diagonal:  bottom(i) -> top(i)
    conn(end+1, :) = [i, top_offset+i]; %#ok<AGROW>
    % Falling diagonal: top(i) -> bottom(i+1)
    conn(end+1, :) = [top_offset+i, i+1]; %#ok<AGROW>
end

num_orig_members = size(conn, 1);

fprintf('=== Original Truss ===\n');
fprintf('  Nodes:   %d  (%d bottom, %d top)\n', n_total, n_bot, n_top);
fprintf('  Members: %d\n', num_orig_members);
fprintf('  Span:    %.0f mm,  Height: %.0f mm\n\n', span*1e3, height*1e3);

%% ==================== SUBDIVIDE ====================

fprintf('\n=== Subdividing (%d sub-elements per member) ===\n', num_sub);

[sub_nodes, sub_conn, member_map] = subdivideTruss(truss_nodes, conn, num_sub);

num_nodes    = size(sub_nodes, 1);
num_elements = size(sub_conn, 1);

fprintf('  Nodes:    %d\n', num_nodes);
fprintf('  Elements: %d\n', num_elements);

%% ==================== PREVIEW ====================

figure('Name', 'Truss Preview', 'Position', [50 500 900 400]);

% -- Left: original truss --
subplot(1,2,1); hold on;
for i = 1:size(conn,1)
    n1 = conn(i,1); n2 = conn(i,2);
    plot(truss_nodes([n1 n2],1)*1e3, truss_nodes([n1 n2],2)*1e3, ...
        'b-', 'LineWidth', 2);
end
plot(truss_nodes(:,1)*1e3, truss_nodes(:,2)*1e3, 'ko', ...
    'MarkerSize', 6, 'MarkerFaceColor', 'k');

% Mark supports
plot(truss_nodes(1,1)*1e3, truss_nodes(1,2)*1e3, ...
    'k^', 'MarkerSize', 14, 'MarkerFaceColor', [0.3 0.3 0.3], 'LineWidth', 2);
plot(truss_nodes(n_bot,1)*1e3, truss_nodes(n_bot,2)*1e3, ...
    'k^', 'MarkerSize', 14, 'MarkerFaceColor', [0.3 0.3 0.3], 'LineWidth', 2);

% Mark load
load_node_orig = ceil(n_bot / 2);   % centre bottom node
quiver(truss_nodes(load_node_orig,1)*1e3, truss_nodes(load_node_orig,2)*1e3, ...
    0, -height*1e3*0.3, 0, 'r', 'LineWidth', 2.5, 'MaxHeadSize', 0.5);

axis equal; grid on;
xlabel('X [mm]'); ylabel('Y [mm]');
title(sprintf('Original Truss (%d members)', num_orig_members));

% -- Right: subdivided mesh --
subplot(1,2,2); hold on;
for i = 1:num_elements
    n1 = sub_conn(i,1); n2 = sub_conn(i,2);
    plot(sub_nodes([n1 n2],1)*1e3, sub_nodes([n1 n2],2)*1e3, ...
        'b-', 'LineWidth', 0.5);
end
plot(sub_nodes(:,1)*1e3, sub_nodes(:,2)*1e3, 'k.', 'MarkerSize', 2);
axis equal; grid on;
xlabel('X [mm]'); ylabel('Y [mm]');
title(sprintf('Subdivided (%d elements)', num_elements));

sgtitle('Truss Definition', 'FontSize', 13, 'FontWeight', 'bold');

%% ==================== FEA SETUP ====================

props = struct();
props.E         = E * ones(num_elements, 1);
props.width     = beam_width * ones(num_elements, 1);
props.thickness = initial_thickness * ones(num_elements, 1);
props.ultimate  = ultimate * ones(num_elements, 1);
props.is_rigid  = false(num_elements, 1);
props.rho       = rho;

% Boundary conditions:
%   Node 1 (bottom-left):  pin  — ux = 0, uy = 0
%   Node n_bot (bottom-right): roller — uy = 0
constraints = [1,     1, 0;
               1,     2, 0;
               n_bot, 2, 0];

geometry = StructureGeometry(sub_nodes, sub_conn, props, constraints);
analysis = FEAAnalysis(geometry);
analysis.verbose = false;

% Load at centre of bottom chord
% The original centre node index is preserved in the subdivided mesh
load_node = load_node_orig;

load_case = LoadCase(geometry);
load_case.addNodalLoad(load_node, 0, Fy_load, 0, 'Centre point load');
load_case.setLoadFactor(1.0);

fprintf('\n=== Loading ===\n');
fprintf('  Fy = %.0f N at node %d  [%.1f, %.1f] mm\n', ...
    Fy_load, load_node, sub_nodes(load_node,:)*1e3);

%% ==================== INITIAL ANALYSIS ====================

fprintf('\n=== Initial Analysis ===\n');

fprintf('  Linear...\n');
lin_results = analysis.runLinearAnalysis(load_case);
lin_stress  = analysis.performStressAnalysis();
fprintf('    Max disp: %.4f mm  |  FOS: %.2f  |  Max stress: %.1f MPa\n', ...
    lin_results.max_displacement*1e3, lin_stress.FOS.min_FOS, ...
    lin_stress.stresses.max_von_mises/1e6);

fprintf('  Nonlinear...\n');
nl_results  = analysis.solveLargeDeflection(load_case, 2);
nl_stress   = analysis.performStressAnalysis();
nl_deformed = analysis.getDeformedCoordinates(1);
fprintf('    Max disp: %.4f mm  |  FOS: %.2f  |  Max stress: %.1f MPa\n', ...
    nl_results.max_displacement*1e3, nl_stress.FOS.min_FOS, ...
    nl_stress.stresses.max_von_mises/1e6);
fprintf('    Mass:     %.4f g\n\n', geometry.getTotalMass()*1e3);

% Store pre-optimisation state
initial_nodes     = geometry.nodes;
initial_thickness = geometry.thickness;
initial_conn      = geometry.connectivity;
initial_mass      = geometry.getTotalMass();

%% ==================== THICKNESS OPTIMISATION ====================

fprintf('=== Thickness Optimisation ===\n');

opt_options = struct();
opt_options.alpha           = 0.3;
opt_options.safety_factor   = 2;
opt_options.max_iterations  = 300;
opt_options.t_min           = 0.0002;       % 0.2 mm  (allow near-zero)
opt_options.t_max           = 0.01;        % 15 mm
opt_options.convergence_tol = 0.0001;
opt_options.move_limit      = 0.15;
opt_options.use_nonlinear   = false;        % linear is faster for trusses
opt_options.verbose         = true;

[optimized_geom, opt_results] = optimizeThickness( ...
    geometry, analysis, load_case, opt_options);

%% ==================== POST-OPTIMISATION ANALYSIS ====================

analysis_opt = FEAAnalysis(optimized_geom);
analysis_opt.verbose = false;
opt_final    = analysis_opt.runLinearAnalysis(load_case);
opt_stress   = analysis_opt.performStressAnalysis();
opt_deformed = analysis_opt.getDeformedCoordinates(1);

fprintf('\n  Post-optimisation:\n');
fprintf('    Max disp: %.4f mm  |  FOS: %.2f  |  Max stress: %.1f MPa\n\n', ...
    opt_final.max_displacement*1e3, opt_stress.FOS.min_FOS, ...
    opt_stress.stresses.max_von_mises/1e6);

%% ==================== PLOTTING ====================

cmap = jet(256);

plotter_init = StructurePlotter(geometry, analysis);
plotter_opt  = StructurePlotter(optimized_geom, analysis_opt);

% ===== FIGURE 1: Initial Deflection + Initial Stress =====
figure('Name', 'Initial Analysis', 'Position', [50 50 1100 475]);

subplot(1,2,1); hold on;
for i = 1:size(initial_conn, 1)
    n1 = initial_conn(i,1); n2 = initial_conn(i,2);
    plotter_init.Plot_Beam( ...
        [initial_nodes(n1,1), initial_nodes(n2,1)], ...
        [initial_nodes(n1,2), initial_nodes(n2,2)], ...
        initial_thickness(i), [0.85 0.85 0.85]);
end
for i = 1:size(initial_conn, 1)
    n1 = initial_conn(i,1); n2 = initial_conn(i,2);
    plotter_init.Plot_Beam( ...
        [nl_deformed(n1,1), nl_deformed(n2,1)], ...
        [nl_deformed(n1,2), nl_deformed(n2,2)], ...
        initial_thickness(i), [0.2 0.4 0.8]);
end
plotter_init.plotSupports();
axis equal; grid on; xlabel('X [m]'); ylabel('Y [m]');
title(sprintf('Initial Deflection (Max: %.3f mm)', nl_results.max_displacement*1e3));

subplot(1,2,2); hold on;
init_vm = nl_stress.stresses.combined;
max_init = max(abs(init_vm(:))) / 1e6;
for i = 1:size(initial_conn, 1)
    n1 = initial_conn(i,1); n2 = initial_conn(i,2);
    s = max(abs(init_vm(i,:))) / 1e6;
    cidx = max(1, min(256, round(s / max(max_init, eps) * 255) + 1));
    plotter_init.Plot_Beam( ...
        [initial_nodes(n1,1), initial_nodes(n2,1)], ...
        [initial_nodes(n1,2), initial_nodes(n2,2)], ...
        initial_thickness(i), cmap(cidx,:));
end
plotter_init.plotSupports();
colormap(gca, cmap); cb = colorbar; cb.Label.String = 'Stress [MPa]';
caxis([0 max_init]);
axis equal; grid on; xlabel('X [m]'); ylabel('Y [m]');
title(sprintf('Initial Stress (Max: %.1f MPa, FOS: %.2f)', ...
    max_init, nl_stress.FOS.min_FOS));

sgtitle('Initial Analysis', 'FontSize', 13, 'FontWeight', 'bold');

% ===== FIGURE 2: Optimised Deflection + Optimised Stress =====
figure('Name', 'Optimised Analysis', 'Position', [50 550 1100 475]);

subplot(1,2,1); hold on;
for i = 1:optimized_geom.num_members
    [n1, n2] = optimized_geom.getMemberNodes(i);
    plotter_opt.Plot_Beam( ...
        [optimized_geom.nodes(n1,1), optimized_geom.nodes(n2,1)], ...
        [optimized_geom.nodes(n1,2), optimized_geom.nodes(n2,2)], ...
        optimized_geom.thickness(i), [0.85 0.85 0.85]);
end
for i = 1:optimized_geom.num_members
    [n1, n2] = optimized_geom.getMemberNodes(i);
    plotter_opt.Plot_Beam( ...
        [opt_deformed(n1,1), opt_deformed(n2,1)], ...
        [opt_deformed(n1,2), opt_deformed(n2,2)], ...
        optimized_geom.thickness(i), [0.8 0.2 0.2]);
end
plotter_opt.plotSupports();
axis equal; grid on; xlabel('X [m]'); ylabel('Y [m]');
title(sprintf('Optimised Deflection (Max: %.3f mm)', opt_final.max_displacement*1e3));

subplot(1,2,2); hold on;
opt_vm = opt_stress.stresses.combined;
max_opt = max(abs(opt_vm(:))) / 1e6;
for i = 1:optimized_geom.num_members
    [n1, n2] = optimized_geom.getMemberNodes(i);
    s = max(abs(opt_vm(i,:))) / 1e6;
    cidx = max(1, min(256, round(s / max(max_opt, eps) * 255) + 1));
    plotter_opt.Plot_Beam( ...
        [optimized_geom.nodes(n1,1), optimized_geom.nodes(n2,1)], ...
        [optimized_geom.nodes(n1,2), optimized_geom.nodes(n2,2)], ...
        optimized_geom.thickness(i), cmap(cidx,:));
end
plotter_opt.plotSupports();
colormap(gca, cmap); cb = colorbar; cb.Label.String = 'Stress [MPa]';
caxis([0 max_opt]);
axis equal; grid on; xlabel('X [m]'); ylabel('Y [m]');
title(sprintf('Optimised Stress (Max: %.1f MPa, FOS: %.2f)', ...
    max_opt, opt_stress.FOS.min_FOS));

sgtitle('Optimised Analysis', 'FontSize', 13, 'FontWeight', 'bold');

% ===== FIGURE 3: Thickness by original member group =====
figure('Name', 'Thickness by Member Group', 'Position', [50 50 900 400]);

subplot(1,2,1); hold on;
% Colour each sub-element by its original member group
group_colors = lines(num_orig_members);
for i = 1:optimized_geom.num_members
    [n1, n2] = optimized_geom.getMemberNodes(i);
    plotter_opt.Plot_Beam( ...
        [optimized_geom.nodes(n1,1), optimized_geom.nodes(n2,1)], ...
        [optimized_geom.nodes(n1,2), optimized_geom.nodes(n2,2)], ...
        optimized_geom.thickness(i), group_colors(member_map(i),:));
end
plotter_opt.plotSupports();
axis equal; grid on; xlabel('X [m]'); ylabel('Y [m]');
title('Thickness by Original Member');

subplot(1,2,2);
bar_data = optimized_geom.thickness * 1e3;
bar(bar_data, 'FaceColor', [0.2 0.4 0.8], 'EdgeColor', 'none');
xlabel('Element'); ylabel('Thickness [mm]');
title('Thickness Distribution');
grid on;
yline(initial_thickness(1)*1e3, 'k--', 'Initial', 'LineWidth', 1.5);
ylim([0 max(bar_data)*1.2]);

sgtitle('Optimised Thickness', 'FontSize', 13, 'FontWeight', 'bold');

%% ==================== EXPORT PROFILE ====================

fprintf('\n=== Exporting Profile ===\n');

profile_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'profiles');
if ~exist(profile_dir, 'dir')
    mkdir(profile_dir);
end

export_opts = struct();
export_opts.num_spline_points = 1000;
export_opts.smoothing         = 0.999999999;
export_opts.scale_to_mm       = 1000;
export_opts.min_thickness_frac = 0.15;    % filter elements < 15% of max
export_opts.export_filename   = fullfile(profile_dir, 'optimised_truss');
export_opts.plot              = true;
export_opts.verbose           = true;

profile = exportProfile(optimized_geom, export_opts);

%% ==================== SUMMARY ====================

fprintf('\n============================\n');
fprintf('  OPTIMISATION SUMMARY\n');
fprintf('============================\n\n');

fprintf('Truss Definition\n');
fprintf('  Span:           %.0f mm,  Height: %.0f mm\n', span*1e3, height*1e3);
fprintf('  Original:       %d members,  %d nodes\n', num_orig_members, n_total);
fprintf('  Subdivided:     %d elements, %d nodes  (%d sub/member)\n', ...
    num_elements, num_nodes, num_sub);
fprintf('  Width:          %.1f mm\n', beam_width*1e3);

fprintf('\nMass\n');
fprintf('  Initial:        %.4f g\n', initial_mass*1e3);
fprintf('  Optimised:      %.4f g\n', opt_results.mass_final*1e3);
fprintf('  Reduction:      %.1f%%\n', opt_results.mass_reduction);

fprintf('\nStress & Safety\n');
fprintf('  Initial FOS:    %.2f\n', nl_stress.FOS.min_FOS);
fprintf('  Optimised FOS:  %.2f\n', opt_stress.FOS.min_FOS);
fprintf('  Max stress:     %.1f MPa\n', opt_stress.stresses.max_von_mises/1e6);
fprintf('  Critical el.:   %d  (original member %d)\n', ...
    opt_stress.stresses.critical_member, ...
    member_map(opt_stress.stresses.critical_member));

fprintf('\nOptimiser\n');
fprintf('  Iterations:     %d\n', opt_results.iterations);
fprintf('  Converged:      %s\n', string(opt_results.converged));
fprintf('  Stress CoV:     %.1f%%\n', opt_results.stress_uniformity*100);

fprintf('\nDeflection\n');
fprintf('  Initial:        %.4f mm\n', nl_results.max_displacement*1e3);
fprintf('  Optimised:      %.4f mm\n', opt_final.max_displacement*1e3);

opt_t = optimized_geom.thickness;
fprintf('\nThickness\n');
fprintf('  Min:            %.2f mm\n', min(opt_t)*1e3);
fprintf('  Max:            %.2f mm\n', max(opt_t)*1e3);
fprintf('  Ratio:          %.1f : 1\n', max(opt_t)/min(opt_t));

% Count effectively removed material
n_thin = sum(opt_t < max(opt_t) * 0.15);
fprintf('  Near-zero:      %d elements (%.0f%%)\n', n_thin, n_thin/num_elements*100);