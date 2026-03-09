classdef StructureGeometry < handle
    % STRUCTUREGEOMETRY  Mutable structure topology, geometry, and material data.
    %
    %   This class stores all geometric, material, and connectivity data for a
    %   2-D frame/truss structure. Properties that change during optimization
    %   (thickness, node positions) are updated through dedicated setter methods
    %   that keep derived quantities (A, I, L, theta0) in sync.
    %
    %   Usage:
    %       props = struct('E', 200e9, 'thickness', 0.003, 'width', 0.01, ...
    %                      'ultimate', 400e6);
    %       constraints = [1 1 0; 1 2 0; 1 3 0];   % fixed node 1
    %       geom = StructureGeometry(nodes, conn, props, constraints);
    %
    %   See also: ElementMatrices, AnalysisState, FEAAnalysis

    % =====================================================================
    %  Topology  (fixed after construction)
    % =====================================================================
    properties (SetAccess = private)
        nodes           % [N x 2]  node coordinates [x, y] (m)
        connectivity    % [M x 2]  element connectivity [start_node, end_node]
    end

    % =====================================================================
    %  Material properties  (per element, fixed after construction)
    % =====================================================================
    properties (SetAccess = private)
        E               % [M x 1]  elastic modulus (Pa)
        ultimate        % [M x 1]  ultimate tensile strength (Pa)
        rho             % [M x 1]  density (kg/m^3)
        width           % [M x 1]  out-of-plane width (m)
        is_rigid        % [M x 1]  logical — treat element as rigid?
    end

    % =====================================================================
    %  Section & geometric properties  (recomputed when thickness or
    %  node positions change)
    % =====================================================================
    properties (SetAccess = private)
        thickness       % [M x 1]  element thickness (m)
        A               % [M x 1]  cross-sectional area  = width .* thickness
        I               % [M x 1]  second moment of area = width .* thickness.^3 / 12
        L               % [M x 1]  element lengths (m)
        theta0          % [M x 1]  initial element orientations (rad)
    end

    % =====================================================================
    %  Boundary conditions
    % =====================================================================
    properties (SetAccess = private)
        constraints     % struct array  .node_id, .dof (1|2|3), .value
    end

    % =====================================================================
    %  Counts
    % =====================================================================
    properties (SetAccess = private)
        num_nodes       % total nodes
        num_members     % total elements
        num_dof_total   % 3 * num_nodes
        num_dof_free    % unconstrained DOFs
        num_constraints % number of constraint equations
    end

    % =====================================================================
    %  DOF & assembly maps  (built once, never change)
    % =====================================================================
    properties (SetAccess = private)
        dof_map         % struct: .free_dof, .fixed_dof, .node_to_dof, .dof_to_node
        assembly_map    % struct: .member_dof [M x 6], .member_nodes [M x 2]
    end

    % =====================================================================
    %  Constructor
    % =====================================================================
    methods
        function obj = StructureGeometry(nodes, connectivity, props, constraints)
            % STRUCTUREGEOMETRY  Create a structure geometry object.
            %
            %   geom = StructureGeometry(nodes, conn, props, constraints)
            %
            %   Inputs
            %   ------
            %   nodes       : [N x 2] node coordinates (m)
            %   connectivity: [M x 2] element connectivity
            %   props       : struct with *required* fields
            %                   .E         — elastic modulus (Pa), scalar or [M x 1]
            %                   .thickness — element thickness (m), scalar or [M x 1]
            %                   .width     — out-of-plane width (m), scalar or [M x 1]
            %                 and *optional* fields
            %                   .ultimate  — UTS (Pa),  default 0 (no strength check)
            %                   .rho       — density (kg/m^3), default 7850 (steel)
            %                   .is_rigid  — logical, default false
            %   constraints : [C x 3] matrix  [node_id, dof, value]
            %                 or struct array with .node_id, .dof, .value

            % --- store topology ---
            obj.nodes        = double(nodes);
            obj.connectivity = int32(connectivity);
            obj.num_nodes    = size(nodes, 1);
            obj.num_members  = size(connectivity, 1);
            obj.num_dof_total = 3 * obj.num_nodes;

            M = obj.num_members;

            % --- required material / section fields ---
            obj.E         = StructureGeometry.expandToVector(props.E, M, 'E');
            obj.thickness = StructureGeometry.expandToVector(props.thickness, M, 'thickness');
            obj.width     = StructureGeometry.expandToVector(props.width, M, 'width');

            % --- optional fields with defaults ---
            if isfield(props, 'ultimate')
                obj.ultimate = StructureGeometry.expandToVector(props.ultimate, M, 'ultimate');
            else
                obj.ultimate = zeros(M, 1);
            end

            if isfield(props, 'rho')
                obj.rho = StructureGeometry.expandToVector(props.rho, M, 'rho');
            else
                obj.rho = 7850 * ones(M, 1);
            end

            if isfield(props, 'is_rigid')
                if isscalar(props.is_rigid)
                    obj.is_rigid = logical(props.is_rigid) & true(M, 1);
                else
                    obj.is_rigid = logical(reshape(props.is_rigid, M, 1));
                end
            else
                obj.is_rigid = false(M, 1);
            end

            % --- derived section properties ---
            obj.A = obj.width .* obj.thickness;
            obj.I = (obj.width .* obj.thickness.^3) / 12;

            % --- constraints ---
            obj.constraints     = StructureGeometry.parseConstraints(constraints);
            obj.num_constraints = numel(obj.constraints);

            % --- geometry ---
            obj.L      = obj.computeLengths();
            obj.theta0 = obj.computeOrientations();

            % --- DOF maps ---
            obj.buildDOFMap();
            obj.num_dof_free = numel(obj.dof_map.free_dof);
            obj.buildAssemblyMap();

            % --- sanity checks ---
            obj.validate();
        end
    end

    % =====================================================================
    %  Mutators  (keep derived quantities in sync)
    % =====================================================================
    methods
        function updateNodePositions(obj, new_nodes)
            % UPDATENODEPOSITIONS  Move nodes and recompute L, theta0.
            %
            %   geom.updateNodePositions(new_nodes)
            %
            %   Input: new_nodes — [N x 2] updated coordinates (m)

            if ~isequal(size(new_nodes), size(obj.nodes))
                error('FEAOPT:badSize', ...
                    'new_nodes must be [%d x 2].', obj.num_nodes);
            end
            obj.nodes  = double(new_nodes);
            obj.L      = obj.computeLengths();
            obj.theta0 = obj.computeOrientations();
        end

        function updateThickness(obj, new_t)
            % UPDATETHICKNESS  Set element thickness and recompute A, I.
            %
            %   geom.updateThickness(new_t)
            %
            %   Input: new_t — [M x 1] thickness vector (m)

            new_t = StructureGeometry.expandToVector(new_t, obj.num_members, 'thickness');
            if any(new_t <= 0)
                error('FEAOPT:nonPositive', 'Thickness must be positive.');
            end
            obj.thickness = new_t;
            obj.A = obj.width .* obj.thickness;
            obj.I = (obj.width .* obj.thickness.^3) / 12;
        end
    end

    % =====================================================================
    %  Queries — element level
    % =====================================================================
    methods
        function [n1, n2] = getMemberNodes(obj, id)
            % GETMEMBERNODES  Node IDs for element id.
            obj.assertMemberID(id);
            n1 = obj.connectivity(id, 1);
            n2 = obj.connectivity(id, 2);
        end

        function dofs = getMemberDOF(obj, id)
            % GETMEMBERDOF  Six global DOF indices for element id  [6 x 1].
            obj.assertMemberID(id);
            dofs = obj.assembly_map.member_dof(id, :)';
        end
    end

    % =====================================================================
    %  Queries — node level
    % =====================================================================
    methods
        function xy = getNodeCoords(obj, id)
            % GETNODECOORDS  [1 x 2] coordinates of node id.
            obj.assertNodeID(id);
            xy = obj.nodes(id, :);
        end

        function dofs = getNodeDOF(obj, id)
            % GETNODEDOF  Three global DOF indices for node id  [3 x 1].
            obj.assertNodeID(id);
            dofs = obj.dof_map.node_to_dof(id, :)';
        end

        function tf = isDOFConstrained(obj, dof_id)
            % ISDOFCONSTRAINED  True if dof_id is fixed.
            tf = ismember(dof_id, obj.dof_map.fixed_dof);
        end

        function [node_id, local_dof] = getDOFNode(obj, dof_id)
            % GETDOFNODE  Owning node and local DOF (1|2|3) for a global DOF.
            if dof_id < 1 || dof_id > obj.num_dof_total
                error('FEAOPT:badDOF', 'DOF %d out of range [1, %d].', ...
                    dof_id, obj.num_dof_total);
            end
            node_id   = obj.dof_map.dof_to_node(dof_id, 1);
            local_dof = obj.dof_map.dof_to_node(dof_id, 2);
        end
    end

    % =====================================================================
    %  Aggregate / convenience queries
    % =====================================================================
    methods
        function m = getMemberMass(obj, id)
            % GETMEMBERMASS  Mass of a single element (kg).
            obj.assertMemberID(id);
            m = obj.rho(id) * obj.A(id) * obj.L(id);
        end

        function m = getTotalMass(obj)
            % GETTOTALMASS  Sum of all element masses (kg).
            m = sum(obj.rho .* obj.A .* obj.L);
        end

        function info = getStructureInfo(obj)
            % GETSTRUCTUREINFO  Summary struct for display / logging.
            info.num_nodes      = obj.num_nodes;
            info.num_members    = obj.num_members;
            info.num_dof_total  = obj.num_dof_total;
            info.num_dof_free   = obj.num_dof_free;
            info.num_constraints = obj.num_constraints;
            info.total_mass     = obj.getTotalMass();
            info.total_length   = sum(obj.L);
        end
    end

    % =====================================================================
    %  Private helpers
    % =====================================================================
    methods (Access = {?AnalysisState, ?ElementMatrices, ?FEAAnalysis, ?StructurePlotter})
        % Bounds-checking helpers — accessible to the core FEA classes.

        function assertMemberID(obj, id)
            if id < 1 || id > obj.num_members
                error('FEAOPT:badMemberID', ...
                    'Member ID %d out of range [1, %d].', id, obj.num_members);
            end
        end

        function assertNodeID(obj, id)
            if id < 1 || id > obj.num_nodes
                error('FEAOPT:badNodeID', ...
                    'Node ID %d out of range [1, %d].', id, obj.num_nodes);
            end
        end
    end

    methods (Access = private)

        % --- geometry computations ---

        function L = computeLengths(obj)
            n1 = obj.connectivity(:,1);
            n2 = obj.connectivity(:,2);
            dx = obj.nodes(n2,1) - obj.nodes(n1,1);
            dy = obj.nodes(n2,2) - obj.nodes(n1,2);
            L  = sqrt(dx.^2 + dy.^2);
        end

        function th = computeOrientations(obj)
            n1 = obj.connectivity(:,1);
            n2 = obj.connectivity(:,2);
            dx = obj.nodes(n2,1) - obj.nodes(n1,1);
            dy = obj.nodes(n2,2) - obj.nodes(n1,2);
            th = atan2(dy, dx);
        end

        % --- DOF map construction ---

        function buildDOFMap(obj)
            N = obj.num_nodes;
            D = obj.num_dof_total;

            % node -> global DOFs
            node_to_dof = zeros(N, 3);
            for i = 1:N
                node_to_dof(i,:) = [3*i-2, 3*i-1, 3*i];
            end

            % global DOF -> (node, local_dof)
            dof_to_node = zeros(D, 2);
            for n = 1:N
                for ld = 1:3
                    gd = 3*(n-1) + ld;
                    dof_to_node(gd,:) = [n, ld];
                end
            end

            % fixed DOFs from constraints
            fixed = zeros(obj.num_constraints, 1);
            for i = 1:obj.num_constraints
                c = obj.constraints(i);
                fixed(i) = 3*(c.node_id - 1) + c.dof;
            end
            fixed = sort(unique(fixed));

            % free DOFs = everything else
            free = setdiff((1:D)', fixed);

            obj.dof_map.node_to_dof = node_to_dof;
            obj.dof_map.dof_to_node = dof_to_node;
            obj.dof_map.fixed_dof   = fixed;
            obj.dof_map.free_dof    = free;
        end

        function buildAssemblyMap(obj)
            M = obj.num_members;
            mdof = zeros(M, 6);
            n1 = obj.connectivity(:,1);
            n2 = obj.connectivity(:,2);
            for i = 1:M
                mdof(i,:) = [3*n1(i)-2, 3*n1(i)-1, 3*n1(i), ...
                             3*n2(i)-2, 3*n2(i)-1, 3*n2(i)];
            end
            obj.assembly_map.member_dof   = mdof;
            obj.assembly_map.member_nodes = obj.connectivity;
        end

        % --- validation ---

        function validate(obj)
            % Topology
            if any(obj.connectivity(:) < 1) || any(obj.connectivity(:) > obj.num_nodes)
                error('FEAOPT:badConnectivity', ...
                    'Connectivity references node IDs outside [1, %d].', obj.num_nodes);
            end
            zeroL = find(obj.L < eps);
            if ~isempty(zeroL)
                error('FEAOPT:zeroLength', ...
                    'Zero-length elements detected: %s', num2str(zeroL'));
            end

            % Constraints
            for i = 1:numel(obj.constraints)
                c = obj.constraints(i);
                if c.node_id < 1 || c.node_id > obj.num_nodes
                    error('FEAOPT:badConstraint', ...
                        'Constraint %d references invalid node %d.', i, c.node_id);
                end
                if c.dof < 1 || c.dof > 3
                    error('FEAOPT:badConstraint', ...
                        'Constraint %d has invalid DOF %d (must be 1, 2, or 3).', i, c.dof);
                end
            end

            if obj.num_dof_free == obj.num_dof_total
                warning('FEAOPT:noConstraints', ...
                    'No constraints applied — rigid-body motion possible.');
            end

            % Material / section
            if any(obj.E <= 0)
                error('FEAOPT:nonPositive', 'Elastic modulus must be positive.');
            end
            if any(obj.A <= 0)
                error('FEAOPT:nonPositive', 'Cross-sectional area must be positive.');
            end
            if any(obj.I <= 0)
                error('FEAOPT:nonPositive', 'Second moment of area must be positive.');
            end
        end

    end

    % =====================================================================
    %  Static utilities
    % =====================================================================
    methods (Static, Access = private)
        function v = expandToVector(val, M, name)
            % Expand a scalar to [M x 1], or validate an existing vector.
            if isscalar(val)
                v = val * ones(M, 1);
            elseif numel(val) == M
                v = reshape(double(val), M, 1);
            else
                error('FEAOPT:badSize', ...
                    '''%s'' must be scalar or length %d.', name, M);
            end
        end

        function cs = parseConstraints(raw)
            % Accept [C x 3] matrix or struct array; return struct array.
            if isnumeric(raw)
                nc = size(raw, 1);
                cs = struct('node_id', num2cell(raw(:,1)), ...
                            'dof',     num2cell(raw(:,2)), ...
                            'value',   num2cell(raw(:,3)));
            elseif isstruct(raw)
                cs = raw;
            else
                error('FEAOPT:badConstraints', ...
                    'Constraints must be a [C x 3] matrix or struct array.');
            end
        end
    end
end