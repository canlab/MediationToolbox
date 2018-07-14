function pts = centroid(varargin)
%CENTROID compute centroid (center of mass) of set of points
%
%   PTS = CENTROID(POINTS, MASS)
%   PTS = CENTROID(PTX, PTY, MASS)
%   Compute center of mass of POINTS, weighte by coefficient MASS. POINTS
%   is a [N*2] array, MASS is [N*1] array, and PTX and PTY are also [N*1]
%   arrays.
%   The result PTS has the same size as MASS.
%
%   PTS = CENTROID(POINTS)  or  PTS = CENTROID(PTX, PTY)
%   Assume equal weights for all points.
%
%   NOTE : for computing polygon centroid, we need to compute centroid of
%   all points inside the polygon, not only polygon vertices.
%
%   TODO : manage also multiple centroid computation in one call.
%   TODO : manage different types of norms
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 07/04/2003.
%

mass = 0;
if nargin==1
    % give only array of points
    var = varargin{1};
    px = var(:,1);
    py = var(:,2);
elseif nargin==2
    % either POINTS+MASS or PX+PY
    var = varargin{1};
    if size(var, 2)>1
        % coefs are POINTS, and MASS
        px = var(:,1);
        py = var(:,2);
        mass = varargin{2};
    else
        % coefs are PX and PY
        px = var;
        py = varargin{2};
    end
elseif nargin==3
    % coefs are PX, PY, and MASS
    px = varargin{1};
    py = varargin{2};
    mass = varargin{3};
end

if mass==0
    mass = ones(length(px), 1);
end
if size(mass,2)>1
    mass=mass';
end


pts(1,1) = mean(px.*mass);
pts(1,2) = mean(py.*mass);

