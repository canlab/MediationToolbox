function varargout = fillPolygon(varargin)
%FILLPOLYGON fill a polygon specified by a list of points
%
%   fillPolygon(COORD) packs coordinates in a single [N*2] array.
%
%   fillPolygon(PX, PY) specify coordinates in separate arrays.
%
%
%   H = fillPolygon(...) also return a handle to the created patch
%
%
%   See also : drawCurve, drawPolygon
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 07/04/2005.
%



color = [0 0 1];
edgeColor = color;

if length(varargin)>0
    var = varargin{1};
    if iscell(var)
        h = zeros(length(var), 1);
        for i=1:length(var)
            h(i) = fillPolygon(var{i}, varargin{2:end});
        end
        return;
    end
end

if length(varargin)==1
    var = varargin{1};

    px = var(:, 1);
    py = var(:, 2);
elseif length(varargin)==2
    var1 = varargin{1};
    var2 = varargin{2};
    
    if size(var1, 1)==size(var2, 1) 
        % same size : parameters are px and py
        px = var1;
        py = var2;
    else
        % different size : second argument is color
        px = var1(:,1);
        py = var1(:,2);
        color = var2;
        edgeColor = color;
    end
elseif length(varargin)>2
    px = varargin{1};
    py = varargin{2};
    color = varargin{3};
    edgeColor = color;
    if length(varargin)>3
        edgeColor = varargin{4};
    end
else
    error ('wrong number of arguments in "fillPolygon"');
end

N = size(px, 1);
if px(N)~=px(1)
    px(size(px, 1)+1, :) = px(1,:);
    py(size(py, 1)+1, :) = py(1,:);
end

h = patch(px, py, color);
set(h, 'edgeColor', edgeColor);


if nargout>0
    varargout{1}=h;
end