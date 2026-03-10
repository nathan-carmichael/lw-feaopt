%% demo_truss_simple.m
%  ShepherdLab LW FEA-OPT — Simple overhanging truss with end moment
%
%  A 4-bay truss pinned at the left and supported at mid-span, with a
%  2-bay overhang. A moment is applied at the free tip to simulate an
%  eccentric attachment or actuator torque.
%
%  See also: subdivideTruss, optimizeThickness, exportProfile

clear; clc; close all;
set(0, 'DefaultAxesFontName', 'Calibri', 'DefaultTextFontName', 'Calibri');

%% ==================== PARAMETERS ====================

% Geometry
bay_width  = 0.060;         % 60 mm per bay
height     = 0.050;         % 50 mm depth
n_bays     = 4;             % total bays (2 supported + 2 overhang)
support_bay = 2;            % right support at bay 2

% Material (steel)
E         = 200e9;
ultimate  = 400e6;
rho       = 7850;

% Cross-section
beam_width        = 0.008;  % 8 mm
initial_thickness = 0.003;  % 3 mm

% Loading
M_tip  = 10;                % 80 Nm moment at free tip
Fy_tip = -3000;              % 500 N downward at free tip

% Subdivision
num_sub = 50;

%% ==================== BUILD TRUSS ====================

% Bottom chord: n_bays+1 nodes at y=0
n_bot = n_bays + 1;
bot_x = (0:n_bays)' * bay_width;
bot_y = zeros(n_bot, 1);

% Top chord: Warren offset at y=height, plus one node directly above
% the left support for the roller constraint
n_top_warren = n_bays;
warren_x = ((0:n_top_warren-1)' + 0.5) * bay_width;
warren_y = height * ones(n_top_warren, 1);

% Extra top node directly above node 1 at (0, height)
top_x = [0; warren_x];
top_y = [height; warren_y];
n_top = numel(top_x);

truss_nodes = [bot_x, bot_y; top_x, top_y];
n_total    = size(truss_nodes, 1);
top_offset = n_bot;

% Node indices
%   bottom:  1 .. n_bot
%   top:     top_offset+1 = node above left support
%            top_offset+2 .. top_offset+n_top = Warren top nodes

% Connectivity
conn = [];

% Bottom chord
for i = 1:n_bot-1
    conn(end+1,:) = [i, i+1]; 
end

% Vertical at left support: bottom node 1 to top node directly above
conn(end+1,:) = [1, top_offset+1]; 

% Top chord: connect the support top node to first Warren node, then chain
for i = 1:n_top-1
    conn(end+1,:) = [top_offset+i, top_offset+i+1]; 
end

% Warren diagonals: connect Warren top nodes (indices 2..n_top in top array)
% to their adjacent bottom nodes
for i = 1:n_top_warren
    top_idx = top_offset + 1 + i;   % skip the support top node
    conn(end+1,:) = [i, top_idx];           %  rising
    conn(end+1,:) = [top_idx, i+1];         %  falling
end

num_orig_members = size(conn, 1);

% Key nodes
pin_node     = 1;                        % left end, bottom
roller_node  = top_offset + 1;           % directly above pin, top
tip_node     = n_bot;                    % free end, bottom

fprintf('=== Simple Overhanging Truss ===\n');
fprintf('  %d bays (%.0f mm each),  height: %.0f mm\n', n_bays, bay_width*1e3, height*1e3);
fprintf('  All joints are rigid (3-DOF frame elements).\n');
fprintf('  Overhang: bays %d-%d\n', support_bay+1, n_bays);
fprintf('  Nodes: %d,  Members: %d\n', n_total, num_orig_members);
fprintf('  Pin: node %d (bottom-left)  |  Fixed: node %d (top-left)\n', ...
    pin_node, roller_node);
fprintf('  Tip: node %d (bottom-right)\n\n', tip_node);

%% ==================== SUBDIVIDE ====================

fprintf('=== Subdividing (%d sub-elements per member) ===\n', num_sub);
[sub_nodes, sub_conn, member_map] = subdivideTruss(truss_nodes, conn, num_sub);
num_nodes    = size(sub_nodes, 1);
num_elements = size(sub_conn, 1);
fprintf('  Nodes: %d,  Elements: %d\n\n', num_nodes, num_elements);

%% ==================== PREVIEW ====================

figure('Name', 'Simple Truss Preview', 'Position', [50 500 800 350]);
hold on;

for i = 1:size(conn,1)
    n1 = conn(i,1); n2 = conn(i,2);
    plot(truss_nodes([n1 n2],1)*1e3, truss_nodes([n1 n2],2)*1e3, ...
        'b-', 'LineWidth', 2);
end
plot(truss_nodes(:,1)*1e3, truss_nodes(:,2)*1e3, 'ko', ...
    'MarkerSize', 6, 'MarkerFaceColor', 'k');

% Supports (both fully fixed)
plot(truss_nodes(pin_node,1)*1e3, truss_nodes(pin_node,2)*1e3, ...
    'k^', 'MarkerSize', 14, 'MarkerFaceColor', [0.3 0.3 0.3], 'LineWidth', 2);
plot(truss_nodes(roller_node,1)*1e3, truss_nodes(roller_node,2)*1e3, ...
    'kv', 'MarkerSize', 14, 'MarkerFaceColor', [0.3 0.3 0.3], 'LineWidth', 2);

% Moment at tip (arc + arrowhead)
tip_xy = truss_nodes(tip_node,:)*1e3;
arc_r  = height*1e3*0.25;
th_arc = linspace(0, 1.5*pi, 40);
plot(tip_xy(1) + arc_r*cos(th_arc), tip_xy(2) + arc_r*sin(th_arc), ...
    'Color', [0.6 0 0.8], 'LineWidth', 2.5);
plot(tip_xy(1) + arc_r*cos(th_arc(end)), tip_xy(2) + arc_r*sin(th_arc(end)), ...
    'v', 'Color', [0.6 0 0.8], 'MarkerSize', 8, 'MarkerFaceColor', [0.6 0 0.8]);

% Force at tip
quiver(tip_xy(1), tip_xy(2), 0, -arc_r*0.8, 0, ...
    'r', 'LineWidth', 2.5, 'MaxHeadSize', 0.5);

% Mark cantilever region (everything right of the fixed supports)
x_first_free = truss_nodes(2,1)*1e3;   % first bottom node after support
x_tip     = truss_nodes(tip_node,1)*1e3;
fill([x_first_free x_tip x_tip x_first_free], ...
     [-height*1e3*0.15 -height*1e3*0.15 height*1e3*1.15 height*1e3*1.15], ...
     'r', 'FaceAlpha', 0.05, 'EdgeColor', 'none');
text(mean([x_first_free x_tip]), height*1e3*1.08, 'Cantilever', ...
    'HorizontalAlignment', 'center', 'Color', [0.7 0 0], 'FontSize', 10);

axis equal; grid on;
xlabel('X [mm]'); ylabel('Y [mm]');
title('Simple Overhanging Truss');
legend({'Members', '', 'Fixed (bottom)', 'Fixed (top)'}, 'Location', 'northwest');

%% ==================== FEA SETUP ====================

props = struct();
props.E         = E * ones(num_elements, 1);
props.width     = beam_width * ones(num_elements, 1);
props.thickness = initial_thickness * ones(num_elements, 1);
props.ultimate  = ultimate * ones(num_elements, 1);
props.is_rigid  = false(num_elements, 1);
props.rho       = rho;

constraints = [pin_node,    1, 0;      % ux = 0
               pin_node,    2, 0;      % uy = 0
               pin_node,    3, 0;      % theta = 0
               roller_node, 1, 0;      % ux = 0
               roller_node, 2, 0;      % uy = 0
               roller_node, 3, 0];     % theta = 0

geometry = StructureGeometry(sub_nodes, sub_conn, props, constraints);
analysis = FEAAnalysis(geometry);
analysis.verbose = false;

load_case = LoadCase(geometry);
load_case.addNodalLoad(tip_node, 0, Fy_tip, M_tip, 'Tip force + moment');
load_case.setLoadFactor(1.0);

fprintf('=== Loading ===\n');
fprintf('  Fixed at bottom-left (node %d): ux=0, uy=0, theta=0\n', pin_node);
fprintf('  Fixed at top-left (node %d): ux=0, uy=0, theta=0\n', roller_node);
fprintf('  Fy = %.0f N,  Mz = %.0f Nm  at node %d (tip)\n\n', ...
    Fy_tip, M_tip, tip_node);

%% ==================== INITIAL ANALYSIS ====================

fprintf('=== Initial Analysis ===\n');
lin_results = analysis.runLinearAnalysis(load_case);
lin_stress  = analysis.performStressAnalysis();
fprintf('  Max disp: %.4f mm  |  FOS: %.2f  |  Max stress: %.1f MPa\n', ...
    lin_results.max_displacement*1e3, lin_stress.FOS.min_FOS, ...
    lin_stress.stresses.max_von_mises/1e6);
fprintf('  Mass:     %.4f g\n\n', geometry.getTotalMass()*1e3);

initial_nodes     = geometry.nodes;
initial_thickness = geometry.thickness;
initial_conn      = geometry.connectivity;
initial_mass      = geometry.getTotalMass();
nl_deformed       = analysis.getDeformedCoordinates(1);
nl_results        = lin_results;
nl_stress         = lin_stress;

%% ==================== OPTIMISE ====================

fprintf('=== Thickness Optimisation ===\n');

opt_options = struct();
opt_options.alpha           = 0.35;
opt_options.safety_factor   = 2;
opt_options.max_iterations  = 300;
opt_options.t_min           = 0.002;
opt_options.t_max           = 0.02;
opt_options.convergence_tol = 0.001;
opt_options.move_limit      = 0.15;
opt_options.use_nonlinear   = false;
opt_options.verbose         = true;

[optimized_geom, opt_results] = optimizeThickness( ...
    geometry, analysis, load_case, opt_options);

%% ==================== POST-OPTIMISATION ====================

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

% ===== FIGURE: Initial + Optimised side by side =====
figure('Name', 'Analysis Comparison', 'Position', [50 50 1200 800]);

% -- Initial deflection --
subplot(2,2,1); hold on;
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
title(sprintf('Initial Deflection (%.3f mm)', nl_results.max_displacement*1e3));

% -- Initial stress --
subplot(2,2,2); hold on;
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
title(sprintf('Initial Stress (%.1f MPa, FOS: %.2f)', max_init, nl_stress.FOS.min_FOS));

% -- Optimised deflection --
subplot(2,2,3); hold on;
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
title(sprintf('Optimised Deflection (%.3f mm)', opt_final.max_displacement*1e3));

% -- Optimised stress --
subplot(2,2,4); hold on;
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
title(sprintf('Optimised Stress (%.1f MPa, FOS: %.2f)', max_opt, opt_stress.FOS.min_FOS));

sgtitle('Simple Overhanging Truss — Before & After', ...
    'FontSize', 13, 'FontWeight', 'bold');

%% ==================== EXPORT ====================

fprintf('\n=== Exporting Profile ===\n');

profile_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'profiles');
if ~exist(profile_dir, 'dir')
    mkdir(profile_dir);
end

export_opts = struct();
export_opts.num_spline_points  = 1000;
export_opts.smoothing          = 0.999999999;
export_opts.scale_to_mm        = 1000;
export_opts.export_filename    = fullfile(profile_dir, 'simple_truss');
export_opts.plot               = true;
export_opts.verbose            = true;

profile = exportProfile(optimized_geom, export_opts);
