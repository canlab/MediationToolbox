function circle = createDirectedCircle(varargin)
%CREATEDIRECTEDCIRCLE create a directed circle
%
%   C = createDirectedCircle(P1, P2, P3) create a circle going through the
%   given points.
%   C is a 1*4 array of the form : [xc yc radius inv].
%   The last parameter is set to 1 if the points are located in clockwise
%   order on the circle.
%
%   C = createDirectedCircle(P1, P2) create the circle whith edge [p1 p2]
%   as diameter
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 12/01/2005.
%

nargs = length(varargin);
x0=0;
y0=0;
r= 0;
if nargs == 2
    % two points
    p1 = varargin{1};
    p2 = varargin{2};
    x0 = (p1(1)+p2(1))/2;
    y0 = (p1(2)+p2(2))/2;
    r = sqrt((p2(1)-p1(1))*(p2(1)-p1(1)) + (p2(2)-p1(2))*(p2(2)-p1(2)))/2;
    d=0;
elseif nargs==3
    % three points
    p1 = varargin{1};
    p2 = varargin{2};
    p3 = varargin{3};

    line1 = medianLine(p1, p2);
    line2 = medianLine(p1, p3);
    point = intersectLines(line1, line2);
    x0 = point(1); y0 = point(2);
    r = sqrt((p1(1)-x0)*(p1(1)-x0) + (p1(2)-y0)*(p1(2)-y0) );
    
    angle = angle3Points(p1, point, p2) + angle3Points(p2, point, p3);
    d = angle>2*pi;
end
    
        
circle = [x0 y0 r d];