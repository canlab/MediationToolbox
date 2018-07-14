function varargout = drawParabola(varargin)
%DRAWPARABOLA draw a parabola on the current axis
%
%   Draw a vertical parabola, defined by its vertex and the distance
%   between vertex and focus (also equal to the distance between vertex 
%   and directrix).
%   Such a parabola admits a vertical axis of symetry.
%
%   The algebraic equation of parabola is :
%      (y-yVertex) = (x-xVertex)^2 / (4*P)
%   A parametric equation of parabola is :
%      x(t) = 2*p*t + xVertex;
%      y(t) = p*t^2 + yVertex;
%
%   Example :
%   drawParabola([50 50 10]);
%   drawParabola([50 50 -5]);
%
%   See Also:
%   drawCircle, drawEllipse
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 02/06/2006.
%

%   HISTORY

if length(varargin)==1
    parabola = varargin{1};
    x0 = parabola(:,1);
    y0 = parabola(:,2);
    p  = parabola(:,3);
elseif length(varargin)>=3
    x0 = varargin{1};
    y0 = varargin{2};
    p  = varargin{3};
else
    error('drawParabola : please specify center x, center y and parameter');
end


t = -100:.1:100;
h = zeros(size(x0));
for i=1:length(x0)
    xt = x0(i) + 2*p(i)*t;
    yt = y0(i) + p(i)*t.^2;

    h(i) = line(xt, yt);
end

if nargout>0
    varargout{1}=h;
end

