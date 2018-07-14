function edge2 = clipEdge(edge, window)
%CLIPEDGE clip an edge with a rectangular window
%
%   edge2 = clipEdge(edge,  window);
%   edge : [x1 y1 x2 y2],
%   window : [xmin xmax ; ymin ymax] or [xmin xmax ymin ymax];
%   return :
%   edge2 = [xc1 yc1 xc2 yc2];
%
%   If clipping is null, return [0 0 0 0];
%
%   if edge is a [nx4] array, return an [nx4] array, corresponding to each
%   clipped edge.
%
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 14/05/2005.
%

%   HISTORY
%   08/01/2007 sort points according to position on edge, not to x coord
%       -> this allows to return edges with same orientation a source, and
%       to keep first or end points at the same position if their are not
%       clipped.

% process data input
if size(window, 1)==1
    window = window([1 3 2 4]);
end

% get limits of window
xmin = window(1);
ymin = window(2);
xmax = window(3);
ymax = window(4);


% convert window limits into lines

lineX0 = [xmin ymin xmax-xmin 0];
lineX1 = [xmin ymax xmax-xmin 0];
lineY0 = [xmin ymin 0 ymax-ymin];
lineY1 = [xmax ymin 0 ymax-ymin];


% compute outcodes of each vertex
p11 = edge(:,1)<xmin; p21 = edge(:,3)<xmin;
p12 = edge(:,1)>xmax; p22 = edge(:,3)>xmax;
p13 = edge(:,2)<ymin; p23 = edge(:,4)<ymin;
p14 = edge(:,2)>ymax; p24 = edge(:,4)>ymax;
out1 = [p11 p12 p13 p14];
out2 = [p21 p22 p23 p24];

% detect edges totally inside window -> no clip.
inside = sum(out1 | out2, 2)==0;

% detect edges totally outside window
outside = sum(out1 & out2, 2)>0;

% select edges not totally outside, and process separately edges totally
% inside window
ind = find(~(inside | outside));


edge2 = zeros(size(edge));
edge2(inside, :) = edge(inside, :);


for i=1:length(ind)
    % current edge
    iedge = edge(ind(i), :);
        
    % compute intersection points with each line of bounding window
    px0 = intersectLineEdge(lineX0, iedge);
    px1 = intersectLineEdge(lineX1, iedge);
    py0 = intersectLineEdge(lineY0, iedge);
    py1 = intersectLineEdge(lineY1, iedge);
    
%   % sort points according to increasing absiss
%   points = sortrows([px0; px1; py0; py1; iedge(1:2); iedge(3:4)], 1);
% 	points = points(isfinite(points(:,1)), :);
     
    % create array of points
    points  = [px0; px1; py0; py1; iedge(1:2); iedge(3:4)];
    
    % and remove infinite points (edges parallel to box edges)
	points  = points(isfinite(points(:,1)), :);
    
    % compute position of remaining points
    pos     = linePosition(points, createLine(iedge(1:2), iedge(3:4)));
    
    % sort points according to their position on the line
    [pos I] = sort(pos);
    points  = points(I, :);

    % select points on the boundary of window
    eps=1e-12;
    ins = find( points(:,1)-xmin>=-eps & points(:,2)-ymin>=-eps & ...
                points(:,1)-xmax<=eps & points(:,2)-ymax<=eps);
            
    if length(ins)>1
        edge2(ind(i), :) = [points(ins(1),:) points(ins(end),:)];
    end        
end


