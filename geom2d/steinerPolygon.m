function nodes = steinerPolygon(points)
%STEINERPOLYGON  create a Steiner polygon from a set of vectors
%
%   NODES = steinerPolygon(POINTS)
%
%   Example
%   n = steinerPolygon([1 0;0 1;1 1]);
%   drawPolygon(n);
%
%   See also
%
%
% ------
% Author: David Legland
% e-mail: david.legland@jouy.inra.fr
% Created: 2006-04-28
% Copyright 2006 INRA - CEPIA Nantes - MIAJ (Jouy-en-Josas).

nodes = [0 0];
for i=1:length(points)
    nodes = [nodes; nodes+repmat(points(i,:), [size(nodes, 1) 1])];
end

K = convhull(nodes(:,1), nodes(:,2));
nodes = nodes(K, :);
    