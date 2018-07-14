function b = onCircle(point, circle)
%ONCIRCLE test if a point is located on a given circle.
%
%   B = onCircle(POINT, CIRCLE) 
%   return true if point is located on the circle
%
%   Example :
%   onCircle([1 0], [0 0 1])
%   returns true, whereas
%   onCircle([1 1], [0 0 1])
%   return false
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 07/04/2004.
%

d = sqrt(sum(power(point - circle(:,1:2), 2), 2));
b = abs(d-circle(:,3))<1e-12;
    