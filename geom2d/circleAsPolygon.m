function varargout = circleAsPolygon(circle, N)
%CIRCLEASPOLYGON convert a circle into a series of points
%
%   P = circleAsPolygon(circle, N);
%   convert circle given as [x0 y0 r] into an array of  [Nx2] double,
%   containing x and y values of points. 
%   The polygon is not closed
%
%   See also :
%   drawCircle, drawPolygon
%
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 06/04/2005.
%

%   HISTORY



t = (0:2*pi/N:2*pi*(1-1/N))';
x = circle(1) + circle(3)*cos(t);
y = circle(2) + circle(3)*sin(t);

if nargout==1
    varargout{1}=[x y];
elseif nargout==2
    varargout{1}=x;
    varargout{2}=y;    
end