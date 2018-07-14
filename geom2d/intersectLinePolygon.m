function pi = intersectLinePolygon(line, points)
%INTERSECTLINEPOLYGON get intersection pts between a line and a polygon
%
%   P = intersectLine(LINE, POLY) returns the intersection points of
%   lines LINE with polygon POLY. LINE is a 1x4 array containing parametric
%   representation of the line ([x0 y0 dx dy], see createLine for details).
%   POLY is a Nx2 array containing coordinate of polygon vertices
%   
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 31/10/2003.
%

%   HISTORY

N = length(points);
edges = [points(1:N-1,:) points(2:N, :)];
edges = [edges; points(N,:) points(1, :)];

pi = zeros(0,2);
for i=1:length(edges)
    p0 = intersectLines(line, createLine(edges(i,1:2), edges(i,3:4)));
    if onEdge(p0, edges(i,:))
        pi(size(pi,1)+1, 1:2) = p0;
    end
end

uniquepi(1,:) = pi(1,:);
for i=2:size(pi, 1)
    isunique = 1;
    for j=1:size(uniquepi, 1)
        if(abs(sum(pi(i,:) - uniquepi(j,:))) < 2*eps)
            isunique = 0;
            break;
        end
    end
    if(isunique)
        uniquepi(end+1,:) = pi(i,:);
    end
end
pi = uniquepi;