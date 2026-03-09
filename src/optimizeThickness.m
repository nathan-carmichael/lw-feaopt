function [geometry, results] = optimizeThickness(geometry, analysis, load_case, options)
% OPTIMIZETHICKNESS  Stress-ratio thickness optimization for 2-D frames.
%
%   Iteratively adjusts element thicknesses so that every (non-rigid)
%   element operates near its allowable stress.  Elements that are
%   overstressed get thicker; elements that are understressed get thinner.
%   The result is a more uniform stress distribution and, typically, a
%   lighter structure.
%
%   The update rule for element i at iteration k is:
%
%       t_new(i) = t_old(i) * ( sigma(i) / sigma_allow )^alpha
%
%   where alpha in (0,1] controls aggressiveness (lower = more cautious).
%
%   Inputs
%   ------
%   geometry  : StructureGeometry — modified IN-PLACE
%   analysis  : FEAAnalysis
%   load_case : LoadCase
%   options   : struct (all fields optional)
%       .alpha              Stress-ratio exponent            (0.4)
%       .max_iterations     Iteration cap                    (100)
%       .convergence_tol    Max relative thickness change    (0.01)
%       .move_limit         Per-iteration change fraction    (0.2)
%       .t_min              Minimum thickness [m]            (0.001)
%       .t_max              Maximum thickness [m]            (0.020)
%       .safety_factor      Divisor on ultimate strength     (1.5)
%       .use_nonlinear      Use large-deflection analysis    (false)
%       .num_load_steps     Load steps for nonlinear         (2)
%       .verbose            Print progress table             (true)
%
%   Outputs
%   -------
%   geometry : the same (mutated) StructureGeometry object
%   results  : struct with fields
%       .converged, .iterations, .thickness_initial, .thickness_final,
%       .thickness_history, .mass_initial, .mass_final, .mass_reduction,
%       .mass_history, .max_stress_history, .convergence_history,
%       .final_stress_results, .stress_uniformity
%
%   Example
%   -------
%       opts.alpha = 0.3;  opts.safety_factor = 2;
%       [geom, res] = optimizeThickness(geom, analysis, load_case, opts);
%
%   See also: FEAAnalysis, StructureGeometry, LoadCase

% =========================================================================
%  Options
% =========================================================================
if nargin < 4, options = struct(); end

opt = applyDefaults(options, struct( ...
    'alpha',            0.4,   ...
    'max_iterations',   100,   ...
    'convergence_tol',  0.01,  ...
    'move_limit',       0.2,   ...
    't_min',            0.001, ...
    't_max',            0.020, ...
    'safety_factor',    1.5,   ...
    'use_nonlinear',    false, ...
    'num_load_steps',   2,     ...
    'verbose',          true));

% =========================================================================
%  Initialisation
% =========================================================================
M = geometry.num_members;
optimizable = ~geometry.is_rigid;
num_opt     = sum(optimizable);

if num_opt == 0
    error('FEAOPT:noOptimizable', 'No non-rigid elements to optimise.');
end

sigma_allow = geometry.ultimate / opt.safety_factor;

t_current = geometry.thickness;

% Pre-allocate history
t_history       = zeros(M, opt.max_iterations);
mass_hist       = zeros(opt.max_iterations, 1);
stress_hist     = zeros(opt.max_iterations, 1);
conv_hist       = zeros(opt.max_iterations, 1);

% Suppress analysis chatter during optimisation
saved_verbose    = analysis.verbose;
analysis.verbose = false;

if opt.verbose
    fprintf('\n=== Thickness Optimisation ===\n');
    fprintf('  Optimisable elements: %d / %d\n', num_opt, M);
    fprintf('  alpha = %.2f,  FOS target = %.2f,  move limit = %.0f%%\n', ...
        opt.alpha, opt.safety_factor, opt.move_limit * 100);
    fprintf('  Thickness bounds: [%.2f, %.2f] mm\n', opt.t_min*1e3, opt.t_max*1e3);
    fprintf('\n  %-5s  %-11s  %-12s  %-10s  %s\n', ...
        'Iter', 'Mass [kg]', 'Max Stress', 'Max dT [%]', 'Status');
    fprintf('  %s\n', repmat('-', 1, 58));
end

% =========================================================================
%  Main loop
% =========================================================================
converged = false;

for iter = 1:opt.max_iterations

    % --- record current thickness ---
    t_history(:, iter) = t_current;

    % --- push thickness into geometry & element matrices ---
    geometry.updateThickness(t_current);
    analysis.element_matrices.resetToInitialConfiguration();

    % --- run analysis ---
    if opt.use_nonlinear && iter > 1
        analysis.solveLargeDeflection(load_case, opt.num_load_steps);
    else
        analysis.runLinearAnalysis(load_case);
    end

    % --- stresses ---
    stress_res   = analysis.performStressAnalysis();
    vm           = stress_res.stresses.von_mises;        % [M x 2]
    max_vm       = max(abs(vm), [], 2);                  % [M x 1]

    % --- stress ratios (non-rigid only) ---
    stress_ratio = zeros(M, 1);
    for i = 1:M
        if optimizable(i)
            if max_vm(i) > eps
                stress_ratio(i) = max_vm(i) / sigma_allow(i);
            else
                stress_ratio(i) = 0.1;      % near-zero stress → shrink
            end
        end
    end

    % --- thickness update ---
    t_new = t_current;
    for i = 1:M
        if ~optimizable(i), continue; end

        t_proposed = t_current(i) * stress_ratio(i)^opt.alpha;

        % move limit
        dt_max     = t_current(i) * opt.move_limit;
        t_proposed = max(t_proposed, t_current(i) - dt_max);
        t_proposed = min(t_proposed, t_current(i) + dt_max);

        % absolute bounds
        t_new(i) = max(opt.t_min, min(opt.t_max, t_proposed));
    end

    % --- convergence metric ---
    rel_change     = abs(t_new - t_current) ./ t_current;
    max_change     = max(rel_change(optimizable));

    % --- bookkeeping ---
    current_mass        = sum(geometry.rho .* geometry.width .* t_current .* geometry.L);
    mass_hist(iter)     = current_mass;
    stress_hist(iter)   = max(max_vm(optimizable));
    conv_hist(iter)     = max_change;

    % --- verbose output ---
    if opt.verbose
        if max_change < opt.convergence_tol
            tag = 'CONVERGED';
        elseif any(stress_ratio(optimizable) > 1)
            tag = sprintf('%d overstressed', sum(stress_ratio(optimizable) > 1));
        else
            tag = 'Optimising';
        end
        fprintf('  %-5d  %-11.4f  %-12.2e  %-10.2f  %s\n', ...
            iter, current_mass, stress_hist(iter), max_change*100, tag);
    end

    % --- update thickness for next iteration (or final state) ---
    t_current = t_new;

    % --- check convergence ---
    if max_change < opt.convergence_tol
        converged = true;
        if opt.verbose
            fprintf('\n  Converged after %d iterations.\n', iter);
        end
        break
    end
end

num_iters = min(iter, opt.max_iterations);

% =========================================================================
%  Final state
% =========================================================================

% Apply final thickness
geometry.updateThickness(t_current);
analysis.element_matrices.resetToInitialConfiguration();

% Final analysis
if opt.use_nonlinear
    analysis.solveLargeDeflection(load_case, opt.num_load_steps);
else
    analysis.runLinearAnalysis(load_case);
end
final_stress = analysis.performStressAnalysis();
final_mass   = sum(geometry.rho .* geometry.width .* t_current .* geometry.L);

% Restore verbosity
analysis.verbose = saved_verbose;

% =========================================================================
%  Package results
% =========================================================================
results.converged           = converged;
results.iterations          = num_iters;
results.thickness_initial   = t_history(:, 1);
results.thickness_final     = t_current;
results.thickness_history   = t_history(:, 1:num_iters);
results.mass_initial        = mass_hist(1);
results.mass_final          = final_mass;
results.mass_reduction      = (mass_hist(1) - final_mass) / mass_hist(1) * 100;
results.mass_history        = mass_hist(1:num_iters);
results.max_stress_history  = stress_hist(1:num_iters);
results.convergence_history = conv_hist(1:num_iters);
results.final_stress_results = final_stress;

% Stress uniformity: coefficient of variation across non-rigid elements
opt_stresses = max_vm(optimizable);
opt_stresses = opt_stresses(isfinite(opt_stresses));
if ~isempty(opt_stresses) && mean(opt_stresses) > 0
    results.stress_uniformity = std(opt_stresses) / mean(opt_stresses);
else
    results.stress_uniformity = Inf;
end

% =========================================================================
%  Summary
% =========================================================================
if opt.verbose
    fprintf('\n=== Optimisation Summary ===\n');
    fprintf('  Initial mass:   %.4f kg\n', results.mass_initial);
    fprintf('  Final mass:     %.4f kg\n', results.mass_final);
    fprintf('  Mass reduction: %.1f%%\n',  results.mass_reduction);
    fprintf('  Min FOS:        %.2f\n',    final_stress.FOS.min_FOS);
    fprintf('  Max stress:     %.1f MPa\n', final_stress.stresses.max_von_mises / 1e6);
    fprintf('  Stress CoV:     %.1f%%\n',  results.stress_uniformity * 100);
    if final_stress.FOS.critical_member > 0
        fprintf('  Critical el.:   %d\n', final_stress.FOS.critical_member);
    end
    if ~converged
        fprintf('\n  WARNING: did not converge within %d iterations.\n', opt.max_iterations);
    end
end

end

% =========================================================================
%  Local utilities
% =========================================================================
function s = applyDefaults(s, defaults)
    flds = fieldnames(defaults);
    for i = 1:numel(flds)
        if ~isfield(s, flds{i})
            s.(flds{i}) = defaults.(flds{i});
        end
    end
end