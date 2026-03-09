classdef FEAAnalysis < handle
    % FEAANALYSIS  Main analysis controller for linear and nonlinear FEA.
    %
    %   Supports:
    %     - Linear elastic analysis            (runLinearAnalysis)
    %     - Newton-Raphson large deflection     (solveLargeDeflection)
    %     - Post-processing: stress & FOS       (performStressAnalysis)
    %
    %   Rigid elements are handled via a stiffness penalty in ElementMatrices
    %   (no Lagrange multipliers).
    %
    %   Usage:
    %       analysis = FEAAnalysis(geometry);
    %       results  = analysis.runLinearAnalysis(load_case);
    %       results  = analysis.solveLargeDeflection(load_case, 5);
    %       stress   = analysis.performStressAnalysis();
    %
    %   See also: StructureGeometry, ElementMatrices, AnalysisState, LoadCase

    % =====================================================================
    %  Core components
    % =====================================================================
    properties (SetAccess = private)
        geometry                % StructureGeometry
        element_matrices        % ElementMatrices
        current_state           % AnalysisState
    end

    % =====================================================================
    %  Deformation tracking
    % =====================================================================
    properties (SetAccess = private)
        cumulative_displacements    % [D x 1]  total displacements from original
        deformed_nodes              % [N x 2]  current deformed node positions
    end

    % =====================================================================
    %  Solver settings
    % =====================================================================
    properties
        verbose                 = true      % print iteration info
        max_iterations          = 1000      % max N-R iterations per load step
        force_tolerance         = 1e-2      % normalised force residual
        energy_tolerance        = 1e-6      % |R · du|
        displacement_tolerance  = 1e-6      % normalised displacement increment
    end

    % =====================================================================
    %  Results / history
    % =====================================================================
    properties (SetAccess = private)
        converged_states        % {1 x S} cell — AnalysisState at each load step
        load_factors            % [1 x S] load factor per step
        iteration_history       % [? x 5] [step, iter, |R|, |du|, energy]
        analysis_time           % cumulative wall-clock time (s)
    end

    % =====================================================================
    %  Post-processing cache
    % =====================================================================
    properties (SetAccess = private)
        member_stresses         % struct from calculateMemberStresses
        stress_results          % struct from performStressAnalysis
    end

    % =====================================================================
    %  Constructor
    % =====================================================================
    methods
        function obj = FEAAnalysis(geometry, options)
            % FEAANALYSIS  Create an analysis object for a given geometry.
            %
            %   analysis = FEAAnalysis(geometry)
            %   analysis = FEAAnalysis(geometry, options)
            %
            %   Options is an optional struct whose field names match the
            %   solver-settings properties (e.g. options.verbose = false).

            if ~isa(geometry, 'StructureGeometry')
                error('FEAOPT:badInput', 'Input must be a StructureGeometry object.');
            end

            obj.geometry         = geometry;
            obj.element_matrices = ElementMatrices(geometry);
            obj.current_state    = AnalysisState(geometry);

            obj.cumulative_displacements = zeros(geometry.num_dof_total, 1);
            obj.deformed_nodes           = geometry.nodes;

            % Apply user overrides
            if nargin > 1 && isstruct(options)
                flds = fieldnames(options);
                for i = 1:numel(flds)
                    if isprop(obj, flds{i})
                        obj.(flds{i}) = options.(flds{i});
                    end
                end
            end

            % Initialise storage
            obj.converged_states = {};
            obj.load_factors     = [];
            obj.iteration_history = [];
            obj.analysis_time    = 0;
        end
    end

    % =====================================================================
    %  Linear analysis
    % =====================================================================
    methods
        function results = runLinearAnalysis(obj, load_case)
            % RUNLINEARANALYSIS  One-shot linear elastic solution.
            %
            %   results = analysis.runLinearAnalysis(load_case)

            obj.assertLoadCase(load_case);
            tic;

            if obj.verbose
                fprintf('\n=== Linear Analysis ===\n');
            end

            obj.resetAnalysis();

            % Assemble and partition
            K     = obj.element_matrices.assembleGlobalStiffnessMatrix(obj.geometry);
            P     = load_case.getScaledLoad();
            free  = obj.geometry.dof_map.free_dof;
            K_ff  = K(free, free);
            P_f   = P(free);

            % Solve
            u_f = K_ff \ P_f;

            % Update state
            obj.current_state.updateDisplacements(u_f);
            obj.current_state.f_external = P;

            % Element results
            obj.current_state.extractMemberDisplacements(obj.element_matrices);
            obj.current_state.computeMemberForces(obj.element_matrices);

            results = obj.packageLinearResults();

            elapsed = toc;
            obj.analysis_time = obj.analysis_time + elapsed;

            if obj.verbose
                fprintf('  Solved in %.3f s  |  Max disp: %.3e m\n', ...
                    elapsed, results.max_displacement);
            end
        end
    end

    % =====================================================================
    %  Large-deflection (Newton-Raphson) analysis
    % =====================================================================
    methods
        function results = solveLargeDeflection(obj, load_case, num_steps)
            % SOLVELARGEDEFLECTION  Incremental Newton-Raphson solution.
            %
            %   results = analysis.solveLargeDeflection(load_case, num_steps)

            obj.assertLoadCase(load_case);
            tic;

            num_rigid = sum(obj.geometry.is_rigid);
            if obj.verbose
                fprintf('\n=== Large Deflection Analysis (Newton-Raphson) ===\n');
                fprintf('  Load steps: %d  |  Max iter: %d  |  Tol: %.1e\n', ...
                    num_steps, obj.max_iterations, obj.force_tolerance);
                if num_rigid > 0
                    fprintf('  Rigid elements: %d (penalty method)\n', num_rigid);
                end
                fprintf('  %s\n', repmat('-', 1, 50));
            end

            obj.resetAnalysis();

            free  = obj.geometry.dof_map.free_dof;
            P_ref = load_case.P_ref;
            disp_hist = zeros(num_steps, 1);

            % ----- load stepping loop -----
            for step = 1:num_steps
                lf = step / num_steps;
                obj.load_factors(step) = lf;

                if obj.verbose
                    fprintf('\n  Step %d/%d  (lambda = %.3f)\n', step, num_steps, lf);
                end

                % External forces for this step
                obj.current_state.f_external = P_ref * lf;
                obj.current_state.f_free_ext = obj.current_state.f_external(free);

                % Newton-Raphson iterations
                converged = obj.newtonRaphsonStep(step);

                if ~converged
                    warning('FEAOPT:notConverged', ...
                        'Step %d did not converge after %d iterations (|R| = %.2e).', ...
                        step, obj.max_iterations, obj.current_state.residual_norm);
                end

                % Track displacements
                obj.cumulative_displacements = obj.current_state.u;

                % Element results
                obj.current_state.extractMemberDisplacements(obj.element_matrices);
                obj.current_state.computeMemberForces(obj.element_matrices);

                disp_hist(step) = obj.current_state.getMaxDisplacement();

                % Store snapshot
                snap = AnalysisState(obj.geometry);
                snap.copyFrom(obj.current_state);
                obj.converged_states{step} = snap;

                if obj.verbose
                    fprintf('    -> Max disp = %.4e m\n', disp_hist(step));
                end
            end

            obj.updateDeformedNodes();

            results = obj.packageNonlinearResults(disp_hist);

            elapsed = toc;
            obj.analysis_time = obj.analysis_time + elapsed;

            if obj.verbose
                fprintf('\n  Analysis complete in %.3f s\n', elapsed);
                fprintf('  Final max disp: %.4e m\n', results.max_displacement);
            end
        end
    end

    % =====================================================================
    %  Stress & factor-of-safety post-processing
    % =====================================================================
    methods
        function results = performStressAnalysis(obj)
            % PERFORMSTRESSANALYSIS  Compute stresses and factor of safety.
            %
            %   results = analysis.performStressAnalysis()
            %
            %   Returns a struct with fields .stresses and .FOS.
            %   Also caches the result in obj.stress_results.

            stresses = obj.calculateMemberStresses();
            fos      = obj.calculateFactorOfSafety(stresses);

            results.stresses = stresses;
            results.FOS      = fos;
            obj.stress_results = results;
        end
    end

    % =====================================================================
    %  Utility / query
    % =====================================================================
    methods
        function resetAnalysis(obj)
            % RESETANALYSIS  Return everything to the undeformed state.

            obj.current_state.reset();
            obj.element_matrices.resetToInitialConfiguration();
            obj.cumulative_displacements = zeros(obj.geometry.num_dof_total, 1);
            obj.deformed_nodes   = obj.geometry.nodes;
            obj.converged_states = {};
            obj.load_factors     = [];
            obj.iteration_history = [];
            obj.member_stresses  = [];
            obj.stress_results   = [];
        end

        function coords = getDeformedCoordinates(obj, scale)
            % GETDEFORMEDCOORDINATES  Node positions including scaled displacements.
            %
            %   coords = analysis.getDeformedCoordinates()        % scale = 1
            %   coords = analysis.getDeformedCoordinates(10)      % exaggerated

            if nargin < 2, scale = 1; end

            coords = obj.geometry.nodes;
            for i = 1:obj.geometry.num_nodes
                dof = obj.geometry.getNodeDOF(i);
                coords(i, 1) = coords(i, 1) + scale * obj.current_state.u(dof(1));
                coords(i, 2) = coords(i, 2) + scale * obj.current_state.u(dof(2));
            end
        end

        function plotConvergence(obj)
            % PLOTCONVERGENCE  Semilogy plots of N-R convergence history.

            if isempty(obj.iteration_history)
                warning('FEAOPT:noHistory', 'No iteration history to plot.');
                return
            end

            figure('Name', 'Convergence History', 'Position', [100 100 1200 400]);
            steps  = unique(obj.iteration_history(:,1));
            colors = lines(numel(steps));
            labels = arrayfun(@(s) sprintf('Step %d', s), steps, 'Uniform', false);

            titles  = {'Force Residual', 'Displacement Increment', 'Energy'};
            cols    = [3, 4, 5];
            tols    = [obj.force_tolerance, obj.displacement_tolerance, obj.energy_tolerance];
            ylabels = {'|R|', '|du|', '|R . du|'};

            for p = 1:3
                subplot(1, 3, p); hold on;
                for k = 1:numel(steps)
                    rows = obj.iteration_history(:,1) == steps(k);
                    semilogy(obj.iteration_history(rows, 2), ...
                             obj.iteration_history(rows, cols(p)), ...
                             '.-', 'Color', colors(k,:), 'LineWidth', 1.5, 'MarkerSize', 8);
                end
                yline(tols(p), 'r--', 'Tolerance');
                xlabel('Iteration'); ylabel(ylabels{p});
                title(titles{p}); grid on;
                if p == 1, legend(labels, 'Location', 'best'); end
            end

            sgtitle('Newton-Raphson Convergence');
        end
    end

    % =====================================================================
    %  Newton-Raphson internals
    % =====================================================================
    methods (Access = private)

        function converged = newtonRaphsonStep(obj, step)
            % One complete set of N-R iterations for a single load step.

            free = obj.geometry.dof_map.free_dof;
            converged = false;

            for iter = 1:obj.max_iterations

                % 1. Update element geometry from current displacements
                obj.updateGeometryFromDisplacements();

                % 2. Internal forces and residual
                obj.current_state.computeInternalForces(obj.element_matrices);
                obj.current_state.computeResidual();

                % 3. Tangent stiffness
                K_t  = obj.element_matrices.assembleGlobalStiffnessMatrix(obj.geometry);
                K_ff = K_t(free, free);

                % 4. Solve for increment
                du = K_ff \ obj.current_state.residual;

                % 5. Check for numerical blow-up
                if any(~isfinite(du))
                    warning('FEAOPT:numericalFailure', ...
                        'NaN/Inf in du at step %d, iteration %d.', step, iter);
                    break
                end

                % 6. Convergence check
                [converged, ~] = obj.checkConvergence(iter, du);

                % 7. Log
                energy = abs(obj.current_state.residual' * du);
                obj.iteration_history(end+1, :) = ...
                    [step, iter, obj.current_state.residual_norm, norm(du), energy];

                if obj.verbose && (iter == 1 || mod(iter, 5) == 0 || converged)
                    fprintf('    Iter %3d:  |R|=%.2e  |du|=%.2e  E=%.2e', ...
                        iter, obj.current_state.residual_norm, norm(du), energy);
                    if converged, fprintf('  [CONVERGED]'); end
                    fprintf('\n');
                end

                if converged, break; end

                % 8. Apply increment for next iteration
                obj.current_state.applyIncrement(du);
            end
        end

        function updateGeometryFromDisplacements(obj)
            % Recompute deformed node positions, then update element
            % lengths and orientations in ElementMatrices.

            geom = obj.geometry;
            u    = obj.current_state.u;

            % Deformed node coordinates
            for n = 1:geom.num_nodes
                dof = geom.getNodeDOF(n);
                obj.deformed_nodes(n, 1) = geom.nodes(n, 1) + u(dof(1));
                obj.deformed_nodes(n, 2) = geom.nodes(n, 2) + u(dof(2));
            end

            % New element geometry
            new_L     = zeros(geom.num_members, 1);
            new_theta = zeros(geom.num_members, 1);

            for m = 1:geom.num_members
                [n1, n2] = geom.getMemberNodes(m);
                dx = obj.deformed_nodes(n2, 1) - obj.deformed_nodes(n1, 1);
                dy = obj.deformed_nodes(n2, 2) - obj.deformed_nodes(n1, 2);
                new_theta(m) = atan2(dy, dx);
                new_L(m)     = sqrt(dx^2 + dy^2);
            end

            obj.element_matrices.updateAllMemberGeometry(new_L, new_theta);
        end

        function [converged, info] = checkConvergence(obj, iter, du)
            % Evaluate force, displacement, and energy convergence criteria.

            R_norm = obj.current_state.residual_norm;

            % Reference values for normalisation
            f_ref = max(norm(obj.current_state.f_free_ext), 1);
            u_ref = max(obj.current_state.getMaxDisplacement(), 1e-6);

            force_err = R_norm     / f_ref;
            disp_err  = norm(du)   / u_ref;
            energy    = abs(obj.current_state.residual' * du);

            force_ok  = force_err < obj.force_tolerance;
            disp_ok   = disp_err  < obj.displacement_tolerance;
            energy_ok = energy    < obj.energy_tolerance;

            % First iteration: accept on force alone (no meaningful du yet)
            if iter == 1
                converged = force_ok;
            else
                converged = force_ok && (disp_ok || energy_ok);
            end

            info.force_error = force_err;
            info.disp_error  = disp_err;
            info.energy      = energy;
        end

        function updateDeformedNodes(obj)
            % Final deformed-coordinate snapshot (called once after stepping).
            u = obj.current_state.u;
            for n = 1:obj.geometry.num_nodes
                dof = obj.geometry.getNodeDOF(n);
                obj.deformed_nodes(n, 1) = obj.geometry.nodes(n, 1) + u(dof(1));
                obj.deformed_nodes(n, 2) = obj.geometry.nodes(n, 2) + u(dof(2));
            end
        end
    end

    % =====================================================================
    %  Stress & FOS internals
    % =====================================================================
    methods (Access = private)

        function stresses = calculateMemberStresses(obj)
            % Compute axial, bending, shear, and von-Mises stresses for
            % every element.  Rigid elements get NaN.

            % Ensure element forces are populated
            if all(obj.current_state.member_forces_local(:) == 0)
                obj.current_state.extractMemberDisplacements(obj.element_matrices);
                obj.current_state.computeMemberForces(obj.element_matrices);
            end

            M = obj.geometry.num_members;

            stresses.axial     = zeros(M, 2);
            stresses.bending   = zeros(M, 2);
            stresses.shear     = zeros(M, 2);
            stresses.von_mises = zeros(M, 2);
            stresses.max_stress = zeros(M, 1);
            stresses.is_rigid  = obj.geometry.is_rigid;

            for i = 1:M
                if obj.geometry.is_rigid(i)
                    stresses.axial(i,:)     = NaN;
                    stresses.bending(i,:)   = NaN;
                    stresses.shear(i,:)     = NaN;
                    stresses.von_mises(i,:) = NaN;
                    stresses.max_stress(i)  = NaN;
                    continue
                end

                Ai = obj.geometry.A(i);
                Ii = obj.geometry.I(i);
                ci = obj.geometry.thickness(i) / 2;   % distance to extreme fibre

                f = obj.current_state.member_forces_local(i, :);
                %   f = [N1, V1, M1, N2, V2, M2]

                stresses.axial(i,:)   = [f(1)/Ai,         f(4)/Ai        ];
                stresses.bending(i,:) = [f(3)*ci/Ii,      f(6)*ci/Ii     ];
                stresses.shear(i,:)   = [1.5*f(2)/Ai,     1.5*f(5)/Ai    ];

                for nd = 1:2
                    sig = stresses.axial(i,nd) + stresses.bending(i,nd);
                    tau = stresses.shear(i,nd);
                    stresses.von_mises(i,nd) = sqrt(sig^2 + 3*tau^2);
                end

                stresses.max_stress(i) = max(stresses.von_mises(i,:));
            end

            % Aggregate maxima (non-rigid only)
            nr = ~obj.geometry.is_rigid;
            stresses = obj.aggregateStressMaxima(stresses, nr);

            % Alias for backward compatibility
            stresses.combined     = stresses.von_mises;
            stresses.max_combined = stresses.max_von_mises;

            obj.member_stresses = stresses;
        end

        function fos = calculateFactorOfSafety(obj, stresses)
            % Factor of safety for every element (NaN for rigid elements).

            M = obj.geometry.num_members;

            fos.stress_based = zeros(M, 2);
            fos.member_FOS   = zeros(M, 1);
            fos.is_rigid     = obj.geometry.is_rigid;

            for i = 1:M
                if obj.geometry.is_rigid(i)
                    fos.stress_based(i,:) = NaN;
                    fos.member_FOS(i)     = NaN;
                    continue
                end

                s_ult = obj.geometry.ultimate(i);

                for nd = 1:2
                    sv = abs(stresses.combined(i, nd));
                    if sv > eps && isfinite(sv)
                        fos.stress_based(i, nd) = s_ult / sv;
                    else
                        fos.stress_based(i, nd) = Inf;
                    end
                end

                fos.member_FOS(i) = min(fos.stress_based(i,:));
            end

            % Global minimum (non-rigid, finite)
            nr = ~obj.geometry.is_rigid;
            vals = fos.member_FOS(nr);
            vals = vals(isfinite(vals));

            if ~isempty(vals)
                [fos.min_FOS, idx] = min(vals);
                nr_ids = find(nr);
                fos.critical_member = nr_ids(idx);
            else
                fos.min_FOS = Inf;
                fos.critical_member = 0;
            end
        end
    end

    % =====================================================================
    %  Results packaging
    % =====================================================================
    methods (Access = private)

        function r = packageLinearResults(obj)
            r.displacements       = obj.current_state.u;
            r.displacements_free  = obj.current_state.u_free;
            r.max_displacement    = obj.current_state.getMaxDisplacement();
            r.external_forces     = obj.current_state.f_external;
            r.member_displ_local  = obj.current_state.member_displ_local;
            r.member_displ_global = obj.current_state.member_displ_global;
            r.nodes_deformed      = obj.getDeformedCoordinates(1);
        end

        function r = packageNonlinearResults(obj, disp_hist)
            r.displacements         = obj.current_state.u;
            r.displacements_free    = obj.current_state.u_free;
            r.max_displacement      = obj.current_state.getMaxDisplacement();
            r.external_forces       = obj.current_state.f_external;
            r.member_displ_local    = obj.current_state.member_displ_local;
            r.member_displ_global   = obj.current_state.member_displ_global;
            r.nodes_deformed        = obj.deformed_nodes;
            r.member_rotations      = obj.element_matrices.theta_current ...
                                      - obj.geometry.theta0;
            r.load_history          = obj.load_factors;
            r.displacement_history  = disp_hist;
            r.num_steps             = numel(obj.load_factors);
            r.state_history         = obj.converged_states;
        end
    end

    % =====================================================================
    %  Static helpers
    % =====================================================================
    methods (Static, Access = private)

        function s = aggregateStressMaxima(s, non_rigid_mask)
            % Compute structure-wide maxima, ignoring rigid/NaN/Inf entries.

            pick = @(vals) vals(isfinite(vals));

            ms = pick(s.max_stress(non_rigid_mask));
            if ~isempty(ms)
                s.max_von_mises = max(ms);
                [~, idx] = max(s.max_stress(non_rigid_mask));
                nr_ids = find(non_rigid_mask);
                s.critical_member = nr_ids(idx);
            else
                s.max_von_mises   = 0;
                s.critical_member = 0;
            end

            ax = pick(abs(s.axial(non_rigid_mask, :)));
            s.max_axial = max(ax(:), [], 'all');
            if isempty(s.max_axial), s.max_axial = 0; end

            bd = pick(abs(s.bending(non_rigid_mask, :)));
            s.max_bending = max(bd(:), [], 'all');
            if isempty(s.max_bending), s.max_bending = 0; end

            sh = pick(abs(s.shear(non_rigid_mask, :)));
            s.max_shear = max(sh(:), [], 'all');
            if isempty(s.max_shear), s.max_shear = 0; end
        end

        function assertLoadCase(lc)
            if ~isa(lc, 'LoadCase')
                error('FEAOPT:badInput', 'Input must be a LoadCase object.');
            end
        end
    end
end