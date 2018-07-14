function varargout = drawPoint(varargin)
%DRAWPOINT draw the point on the axis.
%
%   drawPoint(X, Y) will draw points defined by coordinates X and Y.
%   X and Y are N*1 array, with N being number of points to be drawn.
%   If coordinates of points lie outside the visible area, points are
%   not drawn.
%
%   drawPoint(COORD) packs coordinates in a single [N*2] array.
%
%   drawPoint(..., OPT) will draw each point with given option. OPT is a
%   string compatible with 'plot' model. OPT is a single string, it is not
%   a string array.
%
%
%   H = drawPoint(...) also return a handle to each of the drawn points.
%
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 31/10/2003.
%

%   HISTORY
%   23/02/2004 : add more documentation. Manage different kind of inputs. 
%     Does not draw points outside visible area.


symbol = 'o';

if length(varargin)==1
    var = varargin{1};
    px = var(:, 1);
    py = var(:, 2);
elseif length(varargin)==2
    if ischar(varargin{2})
        var = varargin{1};
        px = var(:,1);
        py = var(:,2);
        symbol = varargin{2};
    else
        px = varargin{1};
        py = varargin{2};
    end
elseif length(varargin)==3
    px = varargin{1};
    py = varargin{2};
    symbol = varargin{3};
    
else
    error ('wrong number of arguments in "drawPoint"');
end


lim = get(gca, 'xlim');
xmin = lim(1);
xmax = lim(2);
lim = get(gca, 'ylim');
ymin = lim(1);
ymax = lim(2);

% check validity for display
ok = px>=xmin;
ok = ok & px<=xmax;
ok = ok & py>=ymin;
ok = ok & py<=ymax;

h = plot(px(ok), py(ok), symbol);

if nargout>0
    varargout{1}=h;
end