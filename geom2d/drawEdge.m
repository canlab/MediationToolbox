function varargout = drawEdge(varargin)
%DRAWEDGE draw the edge given by 2 points
%   
%   usage :
%   drawEdge(x1, y1, x2, y2) draw an edge between the points (x1 y1) and
%   (x2 y2).
%
%   drawEdge([x1 y1 x2 y2]) also works.
%
%   Arguments can be single values or array of size [N*1]. In this case,
%   the function draws multiple edges.
%
%   H = drawEdge(..., OPT), with OPT being a set of pariwise options, can
%   specify color, line width and so on ...
%
%   H = drawEdge(...) return handle(s) to created edges(s)
%
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 31/10/2003.
%

%   HISTORY
%   19/02/2004 : add support for arrays of edges.
%   31/03/2004 : change format of edges to [P1 P2] and variants.
%   28/11/2004 : add support for 3D edges
%   01/08/2005 : add support for drawing options


% default values for parameters
options = {'color', [0 0 1], 'linewidth', 1};
edge = [];
h = [];


% find the number of arguments defining edges
nbVal=0;
for i=1:nargin
    if isnumeric(varargin{i})
        nbVal = nbVal+1;
    else
        % stop at the first non-numeric value
        break;
    end
end

% extract drawing options
options = varargin(nbVal+1:end);


% extract edges characteristics
if nbVal==1
    % all parameters in a single array
    edge = varargin{1};

elseif nbVal==2
    % parameters are two points, or two arrays of points, of size N*2.
    p1 = varargin{1};
    p2 = varargin{2};
    edge = [p1 p2];
    
elseif nbVal==4
    % parameters are 4 parameters of the edge : x1 y1 x2 and y2
    edge = [varargin{1} varargin{2} varargin{3} varargin{4}];
end

% draw the edges
if size(edge, 2)==4
    for i=1:size(edge, 1)
        h(i) = line([edge(i, 1) edge(i, 3)], [edge(i, 2) edge(i, 4)]);
    end
else
    for i=1:size(edge, 1)
        h(i) = line( ...
            [edge(i, 1) edge(i, 4)], ...
            [edge(i, 2) edge(i, 5)], ...
            [edge(i, 3) edge(i, 6)]);
    end
end


if length(options)>0
    for i=1:length(h)
        set(h(i), options{:});
    end
end

if nargout>0
    varargout{1}=h;
end
