classdef LoadCase < handle
    % LOADCASE  Nodal load definition and load-factor management for FEA.
    %
    %   Loads are stored as a reference vector P_ref.  The current applied
    %   load is  P_current = lambda * P_ref,  where lambda is a scalar load
    %   factor that can be incremented for nonlinear load stepping.
    %
    %   Usage:
    %       lc = LoadCase(geometry);
    %       lc.addNodalLoad(5, 0, -1000, 0, 'Tip load');
    %       lc.setLoadFactor(1.0);
    %       P  = lc.getScaledLoad();
    %
    %   See also: StructureGeometry, FEAAnalysis

    % =====================================================================
    %  Load vectors
    % =====================================================================
    properties (SetAccess = private)
        P_ref           % [D x 1]  reference (unscaled) nodal loads
        P_current       % [D x 1]  current loads  = lambda * P_ref
    end

    % =====================================================================
    %  Load factor
    % =====================================================================
    properties (SetAccess = private)
        lambda          % current load factor (scalar)
        lambda_history  % vector of all factors that have been set
    end

    % =====================================================================
    %  Book-keeping
    % =====================================================================
    properties (SetAccess = private)
        load_records    % struct array: .node_id, .fx, .fy, .mz, .description
        geometry        % StructureGeometry reference
        num_nodes
        num_dof_total
        num_dof_free
    end

    % =====================================================================
    %  Constructor
    % =====================================================================
    methods
        function obj = LoadCase(geometry)
            % LOADCASE  Create an empty load case for the given geometry.
            %
            %   lc = LoadCase(geometry)

            if ~isa(geometry, 'StructureGeometry')
                error('FEAOPT:badInput', 'Input must be a StructureGeometry object.');
            end

            obj.geometry      = geometry;
            obj.num_nodes     = geometry.num_nodes;
            obj.num_dof_total = geometry.num_dof_total;
            obj.num_dof_free  = geometry.num_dof_free;

            obj.P_ref     = zeros(obj.num_dof_total, 1);
            obj.P_current = zeros(obj.num_dof_total, 1);

            obj.lambda         = 0;
            obj.lambda_history = [];
            obj.load_records   = struct('node_id', {}, 'fx', {}, 'fy', {}, ...
                                        'mz', {}, 'description', {});
        end
    end

    % =====================================================================
    %  Load definition
    % =====================================================================
    methods
        function addNodalLoad(obj, node_id, fx, fy, mz, description)
            % ADDNODALLOAD  Add forces and/or moment at a node.
            %
            %   lc.addNodalLoad(node_id, fx, fy, mz)
            %   lc.addNodalLoad(node_id, fx, fy, mz, 'description')
            %
            %   Forces are additive: calling this twice on the same node
            %   sums the contributions.

            if nargin < 6, description = sprintf('Load at node %d', node_id); end
            if nargin < 5 || isempty(mz), mz = 0; end
            if nargin < 4 || isempty(fy), fy = 0; end
            if nargin < 3 || isempty(fx), fx = 0; end

            if node_id < 1 || node_id > obj.num_nodes
                error('FEAOPT:badNodeID', ...
                    'Node ID %d out of range [1, %d].', node_id, obj.num_nodes);
            end

            dof = obj.geometry.getNodeDOF(node_id);
            obj.P_ref(dof(1)) = obj.P_ref(dof(1)) + fx;
            obj.P_ref(dof(2)) = obj.P_ref(dof(2)) + fy;
            obj.P_ref(dof(3)) = obj.P_ref(dof(3)) + mz;

            obj.load_records(end+1) = struct('node_id', node_id, ...
                'fx', fx, 'fy', fy, 'mz', mz, 'description', description);

            obj.refreshCurrent();
        end

        function clearLoads(obj)
            % CLEARLOADS  Remove all loads and reset the load factor.

            obj.P_ref(:)       = 0;
            obj.P_current(:)   = 0;
            obj.lambda         = 0;
            obj.lambda_history = [];
            obj.load_records   = struct('node_id', {}, 'fx', {}, 'fy', {}, ...
                                        'mz', {}, 'description', {});
        end
    end

    % =====================================================================
    %  Load factor control
    % =====================================================================
    methods
        function setLoadFactor(obj, lf)
            % SETLOADFACTOR  Set the load factor and update P_current.
            %
            %   lc.setLoadFactor(0.5)

            if lf < 0
                warning('FEAOPT:negativeLambda', ...
                    'Negative load factor (%.3f) applied.', lf);
            end
            obj.lambda = lf;
            obj.lambda_history(end+1) = lf;
            obj.refreshCurrent();
        end

        function incrementLoadFactor(obj, delta)
            % INCREMENTLOADFACTOR  Add delta to the current load factor.
            obj.setLoadFactor(obj.lambda + delta);
        end

        function resetLoadFactor(obj)
            % RESETLOADFACTOR  Set lambda back to zero.
            obj.lambda         = 0;
            obj.lambda_history = [];
            obj.refreshCurrent();
        end
    end

    % =====================================================================
    %  Load retrieval
    % =====================================================================
    methods
        function P = getScaledLoad(obj)
            % GETSCALEDLOAD  Current load vector  (lambda * P_ref).
            P = obj.P_current;
        end

        function Pf = getFreeDOFLoads(obj)
            % GETFREEDOFLOADS  Current loads on free DOFs only.
            Pf = obj.P_current(obj.geometry.dof_map.free_dof);
        end

        function Pf = getReferenceFreeDOFLoads(obj)
            % GETREFERENCEFREEDOFLOADS  Unscaled loads on free DOFs.
            Pf = obj.P_ref(obj.geometry.dof_map.free_dof);
        end

        function applyToState(obj, state)
            % APPLYTOSTATE  Write current loads into an AnalysisState.
            %
            %   lc.applyToState(state)

            if ~isa(state, 'AnalysisState')
                error('FEAOPT:badInput', 'Input must be an AnalysisState object.');
            end
            state.f_external = obj.P_current;
            state.f_free_ext = obj.getFreeDOFLoads();
        end
    end

    % =====================================================================
    %  Load-step generation (convenience for nonlinear analysis)
    % =====================================================================
    methods
        function factors = generateLoadSteps(obj, num_steps, type)
            % GENERATELOADSTEPS  Create a vector of load factors in (0, 1].
            %
            %   factors = lc.generateLoadSteps(10)              % linear
            %   factors = lc.generateLoadSteps(10, 'quadratic') % smaller initial steps
            %
            %   Types:
            %     'linear'    — uniform spacing
            %     'quadratic' — smaller initial increments (good for snap-through)

            if nargin < 3, type = 'linear'; end

            switch lower(type)
                case 'linear'
                    factors = linspace(0, 1, num_steps + 1);
                    factors = factors(2:end);

                case 'quadratic'
                    t = linspace(0, 1, num_steps + 1);
                    factors = t.^2;
                    factors = factors(2:end);

                otherwise
                    error('FEAOPT:badOption', ...
                        'Unknown load-step type ''%s''. Use ''linear'' or ''quadratic''.', type);
            end
        end
    end

    % =====================================================================
    %  Queries
    % =====================================================================
    methods
        function s = getTotalLoad(obj)
            % GETTOTALLOAD  Summary of current applied loads.
            %
            %   s = lc.getTotalLoad()
            %
            %   Returns struct with force_magnitude, moment_magnitude,
            %   force_sum [Fx_total; Fy_total], and moment_sum.

            fx_all = obj.P_current(1:3:end);
            fy_all = obj.P_current(2:3:end);
            mz_all = obj.P_current(3:3:end);

            s.force_magnitude  = norm([fx_all; fy_all]);
            s.moment_magnitude = norm(mz_all);
            s.force_sum        = [sum(fx_all); sum(fy_all)];
            s.moment_sum       = sum(mz_all);
        end

        function [fx, fy, mz] = getNodalLoad(obj, node_id)
            % GETNODALLOAD  Current load components at a single node.

            if node_id < 1 || node_id > obj.num_nodes
                error('FEAOPT:badNodeID', ...
                    'Node ID %d out of range [1, %d].', node_id, obj.num_nodes);
            end
            dof = obj.geometry.getNodeDOF(node_id);
            fx = obj.P_current(dof(1));
            fy = obj.P_current(dof(2));
            mz = obj.P_current(dof(3));
        end

        function tf = hasLoads(obj)
            % HASLOADS  True if any reference load is nonzero.
            tf = any(obj.P_ref ~= 0);
        end
    end

    % =====================================================================
    %  Display
    % =====================================================================
    methods
        function displayLoads(obj)
            % DISPLAYLOADS  Print a formatted summary of all applied loads.

            fprintf('\n=== Load Case Summary ===\n');
            fprintf('  Loads defined: %d\n', numel(obj.load_records));
            fprintf('  Load factor:   %.3f\n\n', obj.lambda);

            if ~isempty(obj.load_records)
                fprintf('  %-4s  %-6s  %12s  %12s  %12s  %s\n', ...
                    '#', 'Node', 'Fx [N]', 'Fy [N]', 'Mz [Nm]', 'Description');
                fprintf('  %s\n', repmat('-', 1, 68));
                for i = 1:numel(obj.load_records)
                    r = obj.load_records(i);
                    fprintf('  %-4d  %-6d  %12.2f  %12.2f  %12.2f  %s\n', ...
                        i, r.node_id, ...
                        r.fx * obj.lambda, ...
                        r.fy * obj.lambda, ...
                        r.mz * obj.lambda, ...
                        r.description);
                end
            end

            s = obj.getTotalLoad();
            fprintf('\n  Totals:  |F| = %.3e N,  |M| = %.3e Nm\n', ...
                s.force_magnitude, s.moment_magnitude);
            fprintf('           SFx = %.3e N,  SFy = %.3e N,  SMz = %.3e Nm\n', ...
                s.force_sum(1), s.force_sum(2), s.moment_sum);
        end
    end

    % =====================================================================
    %  Private
    % =====================================================================
    methods (Access = private)
        function refreshCurrent(obj)
            % Recompute P_current from P_ref and lambda.
            obj.P_current = obj.P_ref * obj.lambda;
        end
    end
end