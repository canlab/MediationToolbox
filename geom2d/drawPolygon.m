function varargout = drawPolygon(varargin)
%DRAWPOLYGON draw a polygon specified by a list of points
%
%   drawPolygon(COORD) packs coordinates in a single [N*2] array.
%
%   drawPolygon(PX, PY) specify coordinates in separate arrays.
%
%
%   H = drawPolygon(...) also return a handle to the list of line objects.
%
%
%   See also : drawCurve
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 05/05/2004.
%

symbol = 'b-';
h = 0;

if isempty(varargin)
    var = varargin{1};    
    if iscell(var)
        for i=1:length(var)
            % check for empty polygons
            if length(var{i})>0
                h(i) = drawPolygon(var{i}, varargin{2:end});
            else
                h(i)=-1;
            end
        end
        
        if nargout>0
            varargout{1}=h;
        end
        
        return;
    end
end


if length(varargin)==1
    var = varargin{1};
    
    % check for empty polygons
    if isempty(var)
        return
    end        
    
    px = var(:, 1);
    py = var(:, 2);
elseif length(varargin)==2
    if ischar(varargin{2})
        var = varargin{1};
        
        % check for empty polygons
        if isempty(var)
            return
        end
        
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
    error ('wrong number of arguments in "drawPolygon"');
end

px(size(px, 1)+1, :) = px(1,:);
py(size(py, 1)+1, :) = py(1,:);

h = plot(px, py, symbol);

if nargout>0
    varargout{1}=h;
end