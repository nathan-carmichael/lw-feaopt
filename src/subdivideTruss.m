function [new_nodes, new_conn, member_map] = subdivideTruss(nodes, connectivity, num_sub)
% SUBDIVIDETRUSS  Split each truss member into sub-elements.
%
%   Takes an original truss mesh and replaces every member with a chain of
%   num_sub smaller elements.  The original joint nodes are preserved; new
%   intermediate nodes are inserted along each member.  The result is a
%   fine mesh suitable for thickness optimisation, where material can
%   redistribute along and between load paths.
%
%   [new_nodes, new_conn, member_map] = subdivideTruss(nodes, conn, num_sub)
%
%   Inputs
%   ------
%   nodes        : [N x 2]  original node coordinates
%   connectivity : [M x 2]  original member connectivity
%   num_sub      : scalar    number of sub-elements per member (>= 1)
%
%   Outputs
%   -------
%   new_nodes   : [N2 x 2]  expanded node array (originals first)
%   new_conn    : [M2 x 2]  new connectivity  (M2 = M * num_sub)
%   member_map  : [M2 x 1]  index of the original member each sub-element
%                            belongs to (useful for grouping / colouring)
%
%   Example
%   -------
%       nodes = [0 0; 1 0; 0.5 0.5];
%       conn  = [1 2; 2 3; 3 1];
%       [n2, c2, map] = subdivideTruss(nodes, conn, 10);
%       % 3 members * 10 = 30 sub-elements, 3 + 3*9 = 30 nodes
%
%   See also: StructureGeometry, optimizeThickness

if num_sub < 1
    error('FEAOPT:badInput', 'num_sub must be >= 1.');
end

num_orig_nodes   = size(nodes, 1);
num_orig_members = size(connectivity, 1);

if num_sub == 1
    % No subdivision — return copies
    new_nodes   = nodes;
    new_conn    = connectivity;
    member_map  = (1:num_orig_members)';
    return
end

% =====================================================================
%  Pre-allocate
% =====================================================================
%   Each member adds (num_sub - 1) intermediate nodes.
%   Total new nodes  = original + M * (num_sub - 1)
%   Total new elems  = M * num_sub

num_new_nodes = num_orig_nodes + num_orig_members * (num_sub - 1);
num_new_elems = num_orig_members * num_sub;

new_nodes   = zeros(num_new_nodes, 2);
new_conn    = zeros(num_new_elems, 2);
member_map  = zeros(num_new_elems, 1);

% Copy original nodes
new_nodes(1:num_orig_nodes, :) = nodes;

next_node = num_orig_nodes + 1;   % index for next intermediate node
next_elem = 1;                     % index for next sub-element

% =====================================================================
%  Subdivide each member
% =====================================================================
for m = 1:num_orig_members
    n1 = connectivity(m, 1);
    n2 = connectivity(m, 2);

    p1 = nodes(n1, :);
    p2 = nodes(n2, :);

    % Parameter values for intermediate points: 1/N, 2/N, ..., (N-1)/N
    % (endpoints are the original nodes)
    fracs = (1:num_sub-1)' / num_sub;

    % Insert intermediate nodes
    int_ids = zeros(num_sub - 1, 1);
    for k = 1:num_sub - 1
        new_nodes(next_node, :) = p1 + fracs(k) * (p2 - p1);
        int_ids(k) = next_node;
        next_node = next_node + 1;
    end

    % Build chain:  n1 -> int(1) -> int(2) -> ... -> int(end) -> n2
    chain = [n1; int_ids; n2];

    for k = 1:num_sub
        new_conn(next_elem, :) = [chain(k), chain(k+1)];
        member_map(next_elem)  = m;
        next_elem = next_elem + 1;
    end
end

end