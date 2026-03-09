classdef AnalysisState < handle
    % ANALYSISSTATE  Mutable state vector container for FEA analysis.
    %
    %   Holds displacement, force, and residual vectors for both the full
    %   DOF system and the reduced (free-DOF) system.  Also stores per-
    %   element end forces and displacements in local and global frames.
    %
    %   The class is designed around the Newton-Raphson workflow:
    %     1.  Set external forces         (f_external, f_free_ext)
    %     2.  Compute internal forces     (computeInternalForces)
    %     3.  Compute residual            (computeResidual)
    %     4.  Solve for increment         (caller solves K * du = R)
    %     5.  Apply increment             (applyIncrement)
    %     6.  Go to 2 until converged
    %
    %   Usage:
    %       state = AnalysisState(geometry);
    %       state.updateDisplacements(u_free);
    %       state.computeInternalForces(elem_matrices);
    %       state.computeResidual();
    %
    %   See also: StructureGeometry, ElementMatrices, FEAAnalysis

    % =====================================================================
    %  Full-system vectors  [D x 1]  (D = 3 * num_nodes)
    % =====================================================================
    properties
        u               % displacement vector (all DOFs)
        f_external      % external force vector
        f_internal      % internal (restoring) force vector
    end

    % =====================================================================
    %  Reduced (free-DOF) vectors  [F x 1]
    % =====================================================================
    properties
        u_free          % free-DOF displacements
        f_free_ext      % external forces on free DOFs
        f_free_int      % internal forces on free DOFs
    end

    % =====================================================================
    %  Newton-Raphson residual
    % =====================================================================
    properties
        residual        % [F x 1]  R = f_free_ext - f_free_int
        residual_norm   % scalar   norm(residual)
        last_delta_u    % [F x 1]  most recent displacement increment
    end

    % =====================================================================
    %  Per-element results  [M x 6]
    % =====================================================================
    properties
        member_forces_local     % end forces   in local  coordinates
        member_forces_global    % end forces   in global coordinates
        member_displ_local      % end displ.   in local  coordinates
        member_displ_global     % end displ.   in global coordinates
    end

    % =====================================================================
    %  Reference geometry & dimensions  (read-only after construction)
    % =====================================================================
    properties (SetAccess = private)
        geometry        % StructureGeometry object
        num_nodes
        num_members
        num_dof_total   % D
        num_dof_free    % F
    end

    % =====================================================================
    %  Constructor
    % =====================================================================
    methods
        function obj = AnalysisState(geometry)
            % ANALYSISSTATE  Create a zero-initialised state for the given geometry.
            %
            %   state = AnalysisState(geometry)

            if ~isa(geometry, 'StructureGeometry')
                error('FEAOPT:badInput', 'Input must be a StructureGeometry object.');
            end

            obj.geometry      = geometry;
            obj.num_nodes     = geometry.num_nodes;
            obj.num_members   = geometry.num_members;
            obj.num_dof_total = geometry.num_dof_total;
            obj.num_dof_free  = geometry.num_dof_free;

            obj.reset();
        end
    end

    % =====================================================================
    %  State management
    % =====================================================================
    methods
        function reset(obj)
            % RESET  Zero every vector and matrix in this state object.

            D = obj.num_dof_total;
            F = obj.num_dof_free;
            M = obj.num_members;

            % Full-system vectors
            obj.u          = zeros(D, 1);
            obj.f_external = zeros(D, 1);
            obj.f_internal = zeros(D, 1);

            % Free-DOF vectors
            obj.u_free     = zeros(F, 1);
            obj.f_free_ext = zeros(F, 1);
            obj.f_free_int = zeros(F, 1);

            % Residual
            obj.residual      = zeros(F, 1);
            obj.residual_norm = 0;
            obj.last_delta_u  = [];

            % Element results
            obj.member_forces_local  = zeros(M, 6);
            obj.member_forces_global = zeros(M, 6);
            obj.member_displ_local   = zeros(M, 6);
            obj.member_displ_global  = zeros(M, 6);
        end

        function copyFrom(obj, src)
            % COPYFROM  Deep-copy every field from another AnalysisState.
            %
            %   state.copyFrom(other_state)

            if ~isa(src, 'AnalysisState')
                error('FEAOPT:badInput', 'Input must be an AnalysisState object.');
            end

            obj.u          = src.u;
            obj.f_external = src.f_external;
            obj.f_internal = src.f_internal;

            obj.u_free     = src.u_free;
            obj.f_free_ext = src.f_free_ext;
            obj.f_free_int = src.f_free_int;

            obj.residual      = src.residual;
            obj.residual_norm = src.residual_norm;
            obj.last_delta_u  = src.last_delta_u;

            obj.member_forces_local  = src.member_forces_local;
            obj.member_forces_global = src.member_forces_global;
            obj.member_displ_local   = src.member_displ_local;
            obj.member_displ_global  = src.member_displ_global;
        end
    end

    % =====================================================================
    %  DOF mapping
    % =====================================================================
    methods
        function updateDisplacements(obj, u_free_new)
            % UPDATEDISPLACEMENTS  Set free-DOF displacements and map to full vector.
            %
            %   state.updateDisplacements(u_free_new)
            %
            %   Constrained DOFs remain zero.

            if numel(u_free_new) ~= obj.num_dof_free
                error('FEAOPT:badSize', ...
                    'Expected length %d, got %d.', obj.num_dof_free, numel(u_free_new));
            end

            obj.u_free = u_free_new(:);
            obj.u(:)   = 0;
            obj.u(obj.geometry.dof_map.free_dof) = obj.u_free;
        end

        function mapFullToFree(obj)
            % MAPFULLTOFREE  Extract free-DOF slices from the full vectors.

            idx = obj.geometry.dof_map.free_dof;
            obj.u_free     = obj.u(idx);
            obj.f_free_ext = obj.f_external(idx);
            obj.f_free_int = obj.f_internal(idx);
        end

        function mapFreeToFull(obj)
            % MAPFREETOFULL  Write free-DOF values back into the full vectors.

            idx = obj.geometry.dof_map.free_dof;
            obj.u(idx)          = obj.u_free;
            obj.f_external(idx) = obj.f_free_ext;
            obj.f_internal(idx) = obj.f_free_int;
        end
    end

    % =====================================================================
    %  Force computation
    % =====================================================================
    methods
        function computeInternalForces(obj, elem_matrices)
            % COMPUTEINTERNALFORCES  Assemble internal (restoring) forces from
            %   the current displacement field using element stiffness matrices.
            %
            %   state.computeInternalForces(elem_matrices)
            %
            %   For each element:  f_global = T' * k_local * T * u_global
            %   Element contributions are assembled into obj.f_internal.
            %
            %   NOTE: Fixed-DOF entries are zeroed after assembly because only
            %   free-DOF forces participate in the Newton-Raphson residual.
            %   This means obj.f_internal does NOT contain reaction forces.
            %   To obtain reactions, compute  R = K * u - f_external  on the
            %   fixed DOFs after convergence.

            geom = obj.geometry;
            obj.f_internal(:) = 0;

            for i = 1:geom.num_members
                dof = geom.getMemberDOF(i);
                u_global = obj.u(dof);

                T_i = elem_matrices.getTransformation(i);
                k_i = elem_matrices.getLocalMatrix(i);

                f_global = T_i' * (k_i * (T_i * u_global));

                obj.f_internal(dof) = obj.f_internal(dof) + f_global;
            end

            % Zero fixed DOFs (see note above)
            obj.f_internal(geom.dof_map.fixed_dof) = 0;

            % Extract free-DOF slice
            obj.f_free_int = obj.f_internal(geom.dof_map.free_dof);
        end

        function computeResidual(obj)
            % COMPUTERESIDUAL  Force imbalance on the free DOFs.
            %
            %   R = f_ext - f_int   (should tend to zero at equilibrium)

            obj.residual      = obj.f_free_ext - obj.f_free_int;
            obj.residual_norm = norm(obj.residual);
        end
    end

    % =====================================================================
    %  Element-level results
    % =====================================================================
    methods
        function extractMemberDisplacements(obj, elem_matrices)
            % EXTRACTMEMBERDISPLACEMENTS  Populate member_displ_local/global
            %   from the current nodal displacement vector.
            %
            %   state.extractMemberDisplacements(elem_matrices)

            geom = obj.geometry;
            for i = 1:obj.num_members
                dof = geom.getMemberDOF(i);
                u_g = obj.u(dof);
                obj.member_displ_global(i,:) = u_g';

                T_i = elem_matrices.getTransformation(i);
                obj.member_displ_local(i,:) = (T_i * u_g)';
            end
        end

        function computeMemberForces(obj, elem_matrices)
            % COMPUTEMEMBERFORCES  Element end forces from current displacements.
            %
            %   state.computeMemberForces(elem_matrices)
            %
            %   Stores results in both local and global frames.

            geom = obj.geometry;
            for i = 1:obj.num_members
                dof = geom.getMemberDOF(i);
                u_g = obj.u(dof);

                T_i = elem_matrices.getTransformation(i);
                k_i = elem_matrices.getLocalMatrix(i);

                u_l = T_i * u_g;
                f_l = k_i * u_l;
                f_g = T_i' * f_l;

                obj.member_forces_local(i,:)  = f_l';
                obj.member_forces_global(i,:) = f_g';
            end
        end
    end

    % =====================================================================
    %  Newton-Raphson increment helpers
    % =====================================================================
    methods
        function applyIncrement(obj, delta_u_free)
            % APPLYINCREMENT  Add a displacement increment to the current state.
            %
            %   state.applyIncrement(delta_u_free)

            obj.last_delta_u = delta_u_free;
            obj.u_free = obj.u_free + delta_u_free;
            obj.u(obj.geometry.dof_map.free_dof) = obj.u_free;
        end

        function scaleState(obj, factor)
            % SCALESTATE  Multiply all displacements by a scalar (for load stepping).
            %
            %   state.scaleState(0.5)

            obj.u      = obj.u      * factor;
            obj.u_free = obj.u_free * factor;
        end
    end

    % =====================================================================
    %  Queries
    % =====================================================================
    methods
        function u_node = getNodalDisplacements(obj, node_id)
            % GETNODALDISPLACEMENTS  [3 x 1] displacement at a node: [ux; uy; theta].
            dof = obj.geometry.getNodeDOF(node_id);   % validation inside getNodeDOF
            u_node = obj.u(dof);
        end

        function theta = getNodalRotations(obj)
            % GETNODALROTATIONS  [N x 1] vector of all nodal rotations.
            theta = obj.u(3:3:obj.num_dof_total);
        end

        function [axial, shear, moment] = getMemberForces(obj, id)
            % GETMEMBERFORCES  Local end forces for an element.
            %
            %   [axial, shear, moment] = state.getMemberForces(id)
            %
            %   Each output is [2 x 1]: values at the start and end nodes.

            obj.geometry.assertMemberID(id);   % reuse geometry validator
            f = obj.member_forces_local(id, :);

            axial  = [f(1); f(4)];
            shear  = [f(2); f(5)];
            moment = [f(3); f(6)];
        end

        function d = getMaxDisplacement(obj)
            % GETMAXDISPLACEMENT  Largest translational displacement magnitude.
            %
            %   Only considers ux and uy (not rotations).

            ux = obj.u(1:3:end);
            uy = obj.u(2:3:end);
            d  = max(sqrt(ux.^2 + uy.^2));
        end
    end

    % =====================================================================
    %  Display
    % =====================================================================
    methods
        function displaySummary(obj)
            % DISPLAYSUMMARY  Print a compact overview of the current state.

            fprintf('\n=== Analysis State Summary ===\n');
            fprintf('Nodes: %d, Elements: %d\n', obj.num_nodes, obj.num_members);
            fprintf('Total DOFs: %d, Free DOFs: %d\n', obj.num_dof_total, obj.num_dof_free);
            fprintf('Max displacement: %.3e m\n', obj.getMaxDisplacement());
            fprintf('Residual norm:    %.3e\n', obj.residual_norm);
        end
    end
end