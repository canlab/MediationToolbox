function varargout = drawCircle(varargin)
%DRAWCIRCLE draw a circle on the current axis
%
%   DRAWCIRCLE(X0, Y0, R) : draw the circle with center (X0,Y0) and the
%   radius R. If X0, Y0 and R are vectors of the same length, draw each
%   circle.
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 31/10/2003.
%

%   HISTORY
%   02/11/2004 : add possibility to draw multiple circles in one call
%   12/01/2005 : allow more than 3 parameters

if length(varargin)==1
    circle = varargin{1};
    x0 = circle(:,1);
    y0 = circle(:,2);
    r  = circle(:,3);
elseif length(varargin)>=3
    x0 = varargin{1};
    y0 = varargin{2};
    r  = varargin{3};
else
    error('drawcircle : please specify center x, center y and radius');
end


t = 0:.01:2*pi;
h = zeros(size(x0));
cot = cos(t);
sit = sin(t);
for i=1:length(x0)
    xt = x0(i) + r(i)*cot;
    yt = y0(i) + r(i)*sit;
    xt = [xt x0(i)+r(i)];
    yt = [yt y0(i)];

    h(i) = line(xt, yt);
end

if nargout>0
    varargout{1}=h;
end