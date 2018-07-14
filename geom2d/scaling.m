function trans = scaling(varargin)
%SCALING return 3*3 matrix of a scale in 2 dimensions
%
%   usage :
%   TRANS = scaling(SX, SY);
%   return the translation corresponding to scaling by SX and SY in the 2
%   main directions.
%   The returned matrix has the form :
%   [SX  0  0]
%   [0  SY  0]
%   [0   0  1]
%
%   TRANS = scaling(SX);
%   Assume SX and SY are equals.
%
%   See also :
%   transformPoint, translation, rotation
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 07/04/2004.


%   HISTORY
%   04/01/2007: rename as scaling

if isempty(varargin)
    sx = 1;
    sy = 1;
elseif length(varargin)==1
    var = varargin{1};
    sx = var(1);
    sy = var(1);
    if length(var)>1
        sy = var(2);
    end
else
    sx = varargin{1};
    sy = varargin{2};
end

trans = [sx 0 0; 0 sy 0; 0 0 1];

