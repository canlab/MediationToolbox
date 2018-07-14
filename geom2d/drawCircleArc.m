function varargout = drawCircleArc(varargin)
%DRAWCIRCLEARC draw a circle arc on the current axis
%
%   drawCircleArc(XC, YC, R, T1, T2) draw ellipse with center (XC, YC),
%   with radius R, between angles T1 and T2 (radians)
%
%   drawCircleArc(PARAM) puts all parameters into one single array.
%
%   H = drawCircleArc(...) returns a handle to the created line object.
%   --------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 12/12/2003.
%

%   HISTORY
%   03/05/2004 : angles are given as radians



if length(varargin)==1
    circle = varargin{1};
    x0 = circle(1);
    y0 = circle(2);
    r  = circle(3);
    t1 = circle(4);
    t2 = circle(5);
elseif length(varargin)==5
    x0 = varargin{1};
    y0 = varargin{2};
    r  = varargin{3};
    t1 = varargin{4};
    t2 = varargin{5};
else
    error('drawcirclearc : please specify center (x and y), radius, start and end angles');
end


h= zeros(length(x0), 1);
for i=1:length(x0)
    t = t1:.01:t2;
    xt = x0(i) + r(i)*cos(t);
    yt = y0(i) + r(i)*sin(t);
    xt = [xt x0(i)+r(i)*cos(t2)];
    yt = [yt y0(i)+r(i)*sin(t2)];

    h(i) = line(xt, yt);
end

if nargout>0
    varargout{1}=h;
end