%% demo_beam_optimization.m
%  ShepherdLab LW FEA-OPT — Beam thickness optimisation demo
%
%  Import a beam curve from file, run FEA, optimise thickness,
%  export the profile, and plot.
%
%  Supported file formats:  .csv, .igs / .iges
%
%  See also: importCSVCurve, importIGSCurve, StructureGeometry,
%            FEAAnalysis, LoadCase, optimizeThickness, exportProfile

clear; clc; close all;
set(0, 'DefaultAxesFontName', 'Calibri', 'DefaultTextFontName', 'Calibri');

%% ==================== PARAMETERS ====================

% Material (Ti-6Al-4V example)
E             = 69e9;           % Young's modulus  [Pa]
ultimate      = 900e6;          % Ultimate tensile strength  [Pa]
rho           = 4820;           % Density  [kg/m^3]

% Cross-section
beam_width        = 0.008;      % Out-of-plane width  [m]
initial_thickness = 0.003;      % Starting thickness  [m]

% Loading
Fx = 0;                         % Horizontal force at free end  [N]
Fy = -333;                      % Vertical force at free end  [N]

% Mesh
num_elements  = 250;

%% ==================== IMPORT GEOMETRY ====================

import_file = 's_test_2.igs';  % <-- change to your filename

import_opts = struct();
import_opts.num_elements = num_elements;
import_opts.scale        = 0.001;       % source file in mm -> metres
import_opts.flip         = false;       % set true to reverse start/end
import_opts.plot         = true;
import_opts.verbose      = true;

[~, ~, ext] = fileparts(import_file);
switch lower(ext)
    case '.csv'
        result = importCSVCurve(import_file, import_opts);
    case {'.igs', '.iges'}
        result = importIGSCurve(import_file, import_opts);
    otherwise
        error('Unsupported format: %s.  Use .csv or .igs/.iges', ext);
end

all_nodes    = result.nodes;
connectivity = result.connectivity;
num_nodes    = result.num_nodes;
num_elements = result.num_elements;

fprintf('\n=== Imported Geometry ===\n');
fprintf('  Arc length:  %.4f mm\n', result.arc_length * 1e3);
fprintf('  Start:       [%.4f, %.4f] mm\n', all_nodes(1,:)*1e3);
fprintf('  End:         [%.4f, %.4f] mm\n', all_nodes(end,:)*1e3);
fprintf('  Start tang:  %.1f deg\n', result.start_angle * 180/pi);
fprintf('  End tang:    %.1f deg\n\n', result.end_angle * 180/pi);

%% ==================== FEA SETUP ====================

fixed_node = 1;
free_node  = num_nodes;

props = struct();
props.E         = E * ones(num_elements, 1);
props.width     = beam_width * ones(num_elements, 1);
props.thickness = initial_thickness * ones(num_elements, 1);
props.ultimate  = ultimate * ones(num_elements, 1);
props.is_rigid  = false(num_elements, 1);
props.rho       = rho;

constraints = [fixed_node, 1, 0;    % ux = 0
               fixed_node, 2, 0;    % uy = 0
               fixed_node, 3, 0;    % theta = 0
               free_node,  1, 0;    % ux = 0  (guided)
               free_node,  3, 0];   % theta = 0

geometry = StructureGeometry(all_nodes, connectivity, props, constraints);
analysis = FEAAnalysis(geometry);
analysis.verbose = false;

fprintf('=== Loading ===\n');
fprintf('  Force: [%.1f, %.1f] N at node %d\n', Fx, Fy, free_node);
fprintf('  Fixed: node %d at [%.4f, %.4f] mm\n\n', ...
    fixed_node, all_nodes(fixed_node,:)*1e3);

load_case = LoadCase(geometry);
load_case.addNodalLoad(free_node, Fx, Fy, 0);
load_case.setLoadFactor(1.0);

%% ==================== INITIAL ANALYSIS ====================

fprintf('=== Initial Analysis ===\n');

fprintf('  Linear...\n');
lin_results = analysis.runLinearAnalysis(load_case);
lin_stress  = analysis.performStressAnalysis();
fprintf('    Max disp:  %.4f mm  |  FOS: %.2f  |  Max stress: %.1f MPa\n', ...
    lin_results.max_displacement*1e3, lin_stress.FOS.min_FOS, ...
    lin_stress.stresses.max_von_mises/1e6);

fprintf('  Nonlinear...\n');
nl_results  = analysis.solveLargeDeflection(load_case, 2);
nl_stress   = analysis.performStressAnalysis();
nl_deformed = analysis.getDeformedCoordinates(1);
fprintf('    Max disp:  %.4f mm  |  FOS: %.2f  |  Max stress: %.1f MPa\n', ...
    nl_results.max_displacement*1e3, nl_stress.FOS.min_FOS, ...
    nl_stress.stresses.max_von_mises/1e6);
fprintf('    Mass:      %.4f g\n\n', geometry.getTotalMass()*1e3);

% Store pre-optimisation state
initial_nodes     = geometry.nodes;
initial_thickness = geometry.thickness;
initial_conn      = geometry.connectivity;
initial_mass      = geometry.getTotalMass();

%% ==================== THICKNESS OPTIMISATION ====================

fprintf('=== Thickness Optimisation ===\n');

opt_options = struct();
opt_options.alpha           = 0.4;
opt_options.safety_factor   = 2;
opt_options.max_iterations  = 200;
opt_options.t_min           = 0.0005;       % 0.5 mm
opt_options.t_max           = 0.020;        % 20 mm
opt_options.convergence_tol = 0.0005;
opt_options.move_limit      = 0.2;
opt_options.use_nonlinear   = true;
opt_options.verbose         = true;

[optimized_geom, opt_results] = optimizeThickness( ...
    geometry, analysis, load_case, opt_options);

%% ==================== POST-OPTIMISATION ANALYSIS ====================

analysis_opt = FEAAnalysis(optimized_geom);
analysis_opt.verbose = false;
opt_final    = analysis_opt.solveLargeDeflection(load_case, 2);
opt_stress   = analysis_opt.performStressAnalysis();
opt_deformed = analysis_opt.getDeformedCoordinates(1);

fprintf('\n  Post-optimisation:\n');
fprintf('    Max disp:  %.4f mm  |  FOS: %.2f  |  Max stress: %.1f MPa\n\n', ...
    opt_final.max_displacement*1e3, opt_stress.FOS.min_FOS, ...
    opt_stress.stresses.max_von_mises/1e6);

%% ==================== PLOTTING ====================

cmap = jet(256);

plotter_init = StructurePlotter(geometry, analysis);
plotter_opt  = StructurePlotter(optimized_geom, analysis_opt);

% ===== FIGURE 1: Initial Deflection + Initial Stress =====
figure('Name', 'Initial Analysis', 'Position', [50 50 1050 475]);

% -- Left: Initial Deflection --
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

% Force arrow on deflection subplot
arrow_len = result.arc_length * 1e3 * 0.08 * 1.5;
force_dir = [Fx, Fy] / max(norm([Fx, Fy]), eps);
quiver(all_nodes(free_node,1), all_nodes(free_node,2), ...
    force_dir(1)*arrow_len*1e-3, force_dir(2)*arrow_len*1e-3, 0, ...
    'r', 'LineWidth', 2.5, 'MaxHeadSize', 0.4);

axis equal; grid on; xlabel('X [m]'); ylabel('Y [m]');
title(sprintf('Initial Deflection (Max: %.3f mm)', nl_results.max_displacement*1e3));

% -- Right: Initial Stress --
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
figure('Name', 'Optimised Analysis', 'Position', [50 550 1050 475]);

% -- Left: Optimised Deflection --
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

% -- Right: Optimised Stress --
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

% ===== FIGURE 3: Thickness Distribution + Stress Bar Chart =====
figure('Name', 'Distributions', 'Position', [50 50 1100 400]);

% -- Left: Thickness --
subplot(1,2,1);
bar_data = optimized_geom.thickness(1:num_elements) * 1e3;
b = bar(1:num_elements, bar_data);
b(1).FaceColor = [0.2 0.4 0.8];
xlabel('Element'); ylabel('Thickness [mm]');
title('Optimised Thickness Distribution');
grid on; ylim([0 max(bar_data(:))*1.2]);
yline(initial_thickness(1)*1e3, 'k--', 'Initial', 'LineWidth', 1.5);

% -- Right: Stress --
subplot(1,2,2); hold on;
bar_stress = opt_stress.stresses.max_stress / 1e6;
bar_stress(isnan(bar_stress)) = 0;
bar(bar_stress, 'FaceColor', [1 0.65 0]);
yline(ultimate/1e6, 'r--', 'LineWidth', 2, 'Label', 'Ultimate');
yline(ultimate/1e6/opt_options.safety_factor, 'b--', 'LineWidth', 1.5, ...
    'Label', sprintf('FOS=%.1f', opt_options.safety_factor));
if opt_stress.stresses.critical_member > 0
    cm = opt_stress.stresses.critical_member;
    plot(cm, bar_stress(cm), 'ko', 'MarkerSize', 8, 'LineWidth', 2);
end
xlabel('Element'); ylabel('Max Von Mises [MPa]');
title('Stress Distribution (Optimised)');
grid on; hold off;

sgtitle('Optimisation Results', 'FontSize', 13, 'FontWeight', 'bold');

%% ==================== EXPORT PROFILE ====================

fprintf('\n=== Exporting Profile ===\n');

% Create profiles output directory
profile_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'profiles');
if ~exist(profile_dir, 'dir')
    mkdir(profile_dir);
end

export_opts = struct();
export_opts.num_spline_points = 1000;
export_opts.scale_to_mm       = 1000;
export_opts.export_filename   = fullfile(profile_dir, 'optimised_beam');
export_opts.plot              = true;
export_opts.verbose           = true;

fprintf('DEBUG: thickness range = [%.4f, %.4f] mm\n', ...
    min(optimized_geom.thickness)*1e3, max(optimized_geom.thickness)*1e3);

profile = exportProfile(optimized_geom, export_opts);

%% ==================== SUMMARY ====================

fprintf('\n============================\n');
fprintf('  OPTIMISATION SUMMARY\n');
fprintf('============================\n\n');

fprintf('Geometry\n');
fprintf('  File:          %s\n', import_file);
fprintf('  Arc length:    %.3f mm\n', result.arc_length*1e3);
fprintf('  Width:         %.1f mm\n', beam_width*1e3);

fprintf('\nMass\n');
fprintf('  Initial:       %.4f g\n', initial_mass*1e3);
fprintf('  Optimised:     %.4f g\n', opt_results.mass_final*1e3);
fprintf('  Reduction:     %.1f%%\n', opt_results.mass_reduction);

fprintf('\nStress & Safety\n');
fprintf('  Initial FOS:   %.2f\n', nl_stress.FOS.min_FOS);
fprintf('  Optimised FOS: %.2f\n', opt_stress.FOS.min_FOS);
fprintf('  Max stress:    %.1f MPa\n', opt_stress.stresses.max_von_mises/1e6);
fprintf('  Ultimate:      %.1f MPa\n', ultimate/1e6);
fprintf('  Critical el.:  %d\n', opt_stress.stresses.critical_member);

fprintf('\nOptimiser\n');
fprintf('  Iterations:    %d\n', opt_results.iterations);
fprintf('  Converged:     %s\n', string(opt_results.converged));
fprintf('  Stress CoV:    %.1f%%\n', opt_results.stress_uniformity*100);

fprintf('\nDeflection\n');
fprintf('  Initial:       %.4f mm\n', nl_results.max_displacement*1e3);
fprintf('  Optimised:     %.4f mm\n', opt_final.max_displacement*1e3);
fprintf('  Change:        %.1f%%\n', ...
    (opt_final.max_displacement - nl_results.max_displacement) / ...
    nl_results.max_displacement * 100);

opt_t = optimized_geom.thickness(1:num_elements);
fprintf('\nThickness\n');
fprintf('  Min:           %.2f mm\n', min(opt_t)*1e3);
fprintf('  Max:           %.2f mm\n', max(opt_t)*1e3);
fprintf('  Ratio:         %.2f : 1\n', max(opt_t)/min(opt_t));

k_init = norm([Fx, Fy]) / nl_results.max_displacement;
k_opt  = norm([Fx, Fy]) / opt_final.max_displacement;
fprintf('\nSpring Rate\n');
fprintf('  Initial:       %.2f N/mm\n', k_init / 1e3);
fprintf('  Optimised:     %.2f N/mm\n', k_opt  / 1e3);