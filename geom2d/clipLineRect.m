function edge = clipLineRect(line, rect)
%CLIPLINERECT clip a line with a polygon
%
%   EDGE = clipLineRect(LINE, RECT);
%   
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 24/07/2006.
%

%   HISTORY

if size(line, 1)==1
    line = repmat(line, size(rect, 1), 1);
elseif size(rect, 1)==1
    rect = repmat(rect, size(line, 1), 1);
elseif size(line, 1) ~= size(rect, 1)
    error('bad sizes for input');
end

edge = zeros(size(line, 1), 4);
for i=1:size(line, 1)
    x = rect(i, 1); y = rect(i, 2); w = rect(i, 3); h = rect(i, 4);
    
	% intersection with axis : x=xmin
	px1 = intersectLineEdge(line(i,:), [x y x+w y]);
	px2 = intersectLineEdge(line(i,:), [x+w y x+w y+h]);
	py1 = intersectLineEdge(line(i,:), [x+w y+h x y+h]);
	py2 = intersectLineEdge(line(i,:), [x y+h x y]);
	
	% sort points along the x coordinate, and  draw a line between
	% the two in the middle
	points = sortrows([px1 ; px2 ; py1 ; py2], 1);
	if points(2,1)>=x-1e-14 && points(2,1)<=x+w+1e-14
        if isfinite(points(3,1))
            edge(i, 1:4) = [points(2,:) points(3,:)];
        else
            edge(i, 1:4) = [points(1,:) points(2,:)];
        end 
    else
        % line outside the rectangle
        edge(i, 1:4) = [0 0 0 0];
	end
end

