classdef ElementMatrices < handle
    % ELEMENTMATRICES  Pre-computed and cached element stiffness/transformation matrices.
    %
    %   Manages local stiffness matrices, coordinate transformation matrices,
    %   and global stiffness matrices for every element. Supports geometry
    %   updates for large-deflection (Newton-Raphson) analysis by tracking
    %   current element lengths and orientations.
    %
    %   Rigid elements are handled via a stiffness penalty method: their
    %   elastic modulus is multiplied by a large factor (default 1e2).
    %
    %   Usage:
    %       em = ElementMatrices(geometry);
    %       K  = em.assembleGlobalStiffnessMatrix(geometry);
    %       em.updateMemberGeometry(3, new_L, new_theta);
    %
    %   See also: StructureGeometry, AnalysisState, FEAAnalysis

    % =====================================================================
    %  References
    % =====================================================================
    properties (SetAccess = private)
        geometry                % StructureGeometry object
    end

    % =====================================================================
    %  Cached matrices  [6 x 6 x M]
    % =====================================================================
    properties (SetAccess = private)
        k_local                 % local stiffness matrices
        T                       % transformation matrices
        k_global                % global stiffness matrices
    end

    % =====================================================================
    %  Current geometric state  (updated during large-deflection analysis)
    % =====================================================================
    properties (SetAccess = private)
        theta_current           % [M x 1]  current element orientations (rad)
        L_current               % [M x 1]  current element lengths (m)
    end

    % =====================================================================
    %  Dirty flags — track which matrices need recomputation
    % =====================================================================
    properties (SetAccess = private)
        dirty_local             % [M x 1]  logical
        dirty_global            % [M x 1]  logical
    end

    % =====================================================================
    %  Settings
    % =====================================================================
    properties (SetAccess = private)
        num_members             % convenience copy of geometry.num_members
        rigid_stiffness_factor  % E multiplier for rigid elements (default 1e2)
    end

    % =====================================================================
    %  Constructor
    % =====================================================================
    methods
        function obj = ElementMatrices(geometry)
            % ELEMENTMATRICES  Construct from a StructureGeometry object.
            %
            %   em = ElementMatrices(geometry)

            if ~isa(geometry, 'StructureGeometry')
                error('FEAOPT:badInput', 'Input must be a StructureGeometry object.');
            end

            obj.geometry    = geometry;
            obj.num_members = geometry.num_members;

            % Default penalty factor for rigid elements
            obj.rigid_stiffness_factor = 1e2;

            % Pre-allocate
            M = obj.num_members;
            obj.k_local  = zeros(6, 6, M);
            obj.T        = zeros(6, 6, M);
            obj.k_global = zeros(6, 6, M);

            % Initialise geometric state to undeformed configuration
            obj.theta_current = geometry.theta0;
            obj.L_current     = geometry.L;

            % Clean flags (everything will be computed now)
            obj.dirty_local  = false(M, 1);
            obj.dirty_global = false(M, 1);

            % Build all matrices from scratch
            obj.rebuildAll();
        end
    end

    % =====================================================================
    %  Full rebuild
    % =====================================================================
    methods
        function rebuildAll(obj)
            % REBUILDALL  Recompute every matrix from current L / theta.
            for i = 1:obj.num_members
                obj.k_local(:,:,i)  = obj.buildLocalMatrix(i);
                obj.T(:,:,i)        = ElementMatrices.buildTransformationMatrix(obj.theta_current(i));
                obj.k_global(:,:,i) = obj.T(:,:,i)' * obj.k_local(:,:,i) * obj.T(:,:,i);
            end
            obj.dirty_local(:)  = false;
            obj.dirty_global(:) = false;
        end

        function resetToInitialConfiguration(obj)
            % RESETTOINITIALCONFIGURATION  Restore original geometry and matrices.
            obj.theta_current = obj.geometry.theta0;
            obj.L_current     = obj.geometry.L;
            obj.rebuildAll();
        end
    end

    % =====================================================================
    %  Geometry updates  (called during Newton-Raphson iterations)
    % =====================================================================
    methods
        function updateMemberGeometry(obj, id, new_L, new_theta)
            % UPDATEMEMBERGEOMETRY  Update length and/or orientation for one element.
            %
            %   em.updateMemberGeometry(id, new_L, new_theta)
            %
            %   Directly recomputes affected matrices rather than deferring,
            %   because Newton-Raphson iterations need them immediately.

            obj.assertID(id);

            L_changed     = abs(new_L     - obj.L_current(id))     / obj.L_current(id) > 1e-8;
            theta_changed = abs(new_theta  - obj.theta_current(id)) > 1e-10;

            if ~L_changed && ~theta_changed
                return
            end

            obj.L_current(id)     = new_L;
            obj.theta_current(id) = new_theta;

            if L_changed
                obj.k_local(:,:,id) = obj.buildLocalMatrix(id);
            end

            if theta_changed
                obj.T(:,:,id) = ElementMatrices.buildTransformationMatrix(new_theta);
            end

            % Global matrix always needs updating if either changed
            Tk = obj.T(:,:,id);
            kl = obj.k_local(:,:,id);
            Kg = Tk' * kl * Tk;
            Kg(abs(Kg) < 1e-12) = 0;          % clean numerical dust
            obj.k_global(:,:,id) = Kg;
        end

        function updateAllMemberGeometry(obj, lengths, thetas)
            % UPDATEALLMEMBERGEOMETRY  Batch update for every element.
            %
            %   em.updateAllMemberGeometry(lengths, thetas)
            %
            %   Inputs:  lengths — [M x 1],  thetas — [M x 1]

            M = obj.num_members;
            if numel(lengths) ~= M || numel(thetas) ~= M
                error('FEAOPT:badSize', 'Input vectors must have length %d.', M);
            end
            for i = 1:M
                obj.updateMemberGeometry(i, lengths(i), thetas(i));
            end
        end
    end

    % =====================================================================
    %  Matrix retrieval  (lazy recomputation via dirty flags)
    % =====================================================================
    methods
        function k = getLocalMatrix(obj, id)
            % GETLOCALMATRIX  Return [6 x 6] local stiffness for element id.
            obj.assertID(id);
            if obj.dirty_local(id)
                obj.k_local(:,:,id)  = obj.buildLocalMatrix(id);
                obj.dirty_local(id)  = false;
            end
            k = obj.k_local(:,:,id);
        end

        function K = getGlobalMatrix(obj, id)
            % GETGLOBALMATRIX  Return [6 x 6] global stiffness for element id.
            obj.assertID(id);
            if obj.dirty_local(id) || obj.dirty_global(id)
                if obj.dirty_local(id)
                    obj.k_local(:,:,id) = obj.buildLocalMatrix(id);
                    obj.dirty_local(id) = false;
                end
                Tk = obj.T(:,:,id);
                kl = obj.k_local(:,:,id);
                Kg = Tk' * kl * Tk;
                Kg(abs(Kg) < 1e-12) = 0;
                obj.k_global(:,:,id) = Kg;
                obj.dirty_global(id) = false;
            end
            K = obj.k_global(:,:,id);
        end

        function Tmat = getTransformation(obj, id)
            % GETTRANSFORMATION  Return [6 x 6] transformation matrix for element id.
            obj.assertID(id);
            Tmat = obj.T(:,:,id);
        end

        function theta = getCurrentAngle(obj, id)
            % GETCURRENTANGLE  Current orientation (rad) for element id.
            obj.assertID(id);
            theta = obj.theta_current(id);
        end

        function Lval = getCurrentLength(obj, id)
            % GETCURRENTLENGTH  Current length (m) for element id.
            obj.assertID(id);
            Lval = obj.L_current(id);
        end
    end

    % =====================================================================
    %  Global assembly
    % =====================================================================
    methods
        function K = assembleGlobalStiffnessMatrix(obj, geometry)
            % ASSEMBLEGLOBALSTIFFNESSMATRIX  Assemble full sparse [D x D] matrix.
            %
            %   K = em.assembleGlobalStiffnessMatrix(geometry)

            % Flush any pending lazy updates
            obj.flushDirty();

            D = geometry.num_dof_total;
            K = sparse(D, D);

            for i = 1:obj.num_members
                dof = geometry.getMemberDOF(i);
                K(dof, dof) = K(dof, dof) + obj.k_global(:,:,i);
            end
        end
    end

    % =====================================================================
    %  Rigid stiffness factor
    % =====================================================================
    methods
        function setRigidStiffnessFactor(obj, factor)
            % SETRIGIDSTIFFNESSFACTOR  Change the E multiplier for rigid elements.
            %
            %   em.setRigidStiffnessFactor(1e4)

            if factor <= 1
                warning('FEAOPT:lowRigidFactor', ...
                    'Rigid stiffness factor should be >> 1 for effective rigid behaviour.');
            end

            obj.rigid_stiffness_factor = factor;

            % Mark rigid elements dirty so they get rebuilt on next access
            rigid_ids = find(obj.geometry.is_rigid);
            obj.dirty_local(rigid_ids)  = true;
            obj.dirty_global(rigid_ids) = true;
        end
    end

    % =====================================================================
    %  Diagnostics
    % =====================================================================
    methods
        function info = getUpdateStatus(obj)
            % GETUPDATESTATUS  How many matrices are flagged dirty.
            info.num_dirty_local  = sum(obj.dirty_local);
            info.num_dirty_global = sum(obj.dirty_global);
            info.members_dirty_local  = find(obj.dirty_local);
            info.members_dirty_global = find(obj.dirty_global);
            info.pct_dirty = 100 * (info.num_dirty_local + info.num_dirty_global) ...
                             / (2 * obj.num_members);
        end
    end

    % =====================================================================
    %  Private — matrix builders
    % =====================================================================
    methods (Access = private)

        function k = buildLocalMatrix(obj, id)
            % Standard 2-D frame element in local coordinates.
            %
            %   DOF order: [u1, v1, theta1, u2, v2, theta2]
            %
            %   Rigid elements get E multiplied by rigid_stiffness_factor.

            Ev = obj.geometry.E(id);
            Av = obj.geometry.A(id);
            Iv = obj.geometry.I(id);
            Lv = obj.L_current(id);

            if obj.geometry.is_rigid(id)
                Ev = Ev * obj.rigid_stiffness_factor;
            end

            ea_L = Ev * Av / Lv;              % axial stiffness  EA/L
            ei_L3 = Ev * Iv / Lv^3;           % flexural coefficient  EI/L^3

            k = [  ea_L,        0,            0,           -ea_L,        0,            0;
                   0,           12*ei_L3,     6*ei_L3*Lv,   0,          -12*ei_L3,     6*ei_L3*Lv;
                   0,           6*ei_L3*Lv,   4*ei_L3*Lv^2, 0,          -6*ei_L3*Lv,   2*ei_L3*Lv^2;
                  -ea_L,        0,            0,            ea_L,        0,            0;
                   0,          -12*ei_L3,    -6*ei_L3*Lv,   0,           12*ei_L3,    -6*ei_L3*Lv;
                   0,           6*ei_L3*Lv,   2*ei_L3*Lv^2, 0,          -6*ei_L3*Lv,   4*ei_L3*Lv^2];
        end

        function flushDirty(obj)
            % Recompute any matrices that are still flagged dirty.
            for i = find(obj.dirty_local)'
                obj.k_local(:,:,i) = obj.buildLocalMatrix(i);
                obj.dirty_local(i) = false;
            end
            for i = find(obj.dirty_global | obj.dirty_local)'
                Tk = obj.T(:,:,i);
                kl = obj.k_local(:,:,i);
                Kg = Tk' * kl * Tk;
                Kg(abs(Kg) < 1e-12) = 0;
                obj.k_global(:,:,i) = Kg;
                obj.dirty_global(i) = false;
            end
        end

        function assertID(obj, id)
            if id < 1 || id > obj.num_members
                error('FEAOPT:badMemberID', ...
                    'Element ID %d out of range [1, %d].', id, obj.num_members);
            end
        end
    end

    % =====================================================================
    %  Static — transformation matrix (pure function of angle)
    % =====================================================================
    methods (Static, Access = private)
        function Tmat = buildTransformationMatrix(theta)
            % 2-D frame rotation: local <-> global for one element.
            c = cos(theta);
            s = sin(theta);
            Tmat = [ c  s  0  0  0  0;
                    -s  c  0  0  0  0;
                     0  0  1  0  0  0;
                     0  0  0  c  s  0;
                     0  0  0 -s  c  0;
                     0  0  0  0  0  1];
        end
    end
end