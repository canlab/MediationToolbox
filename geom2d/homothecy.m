function trans = homothecy(point, ratio)
%HOMOTHECY create a homothecy as an affine transform
%
%   usage :
%   TRANS = HOMOTHECY(POINT, RATIO)
%
%   See also :
%   transformPoint, translation
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 20/01/2005.
%

x0 = point(:,1);
y0 = point(:,2);

m00 = ratio;
m01 = 0;
m02 = x0*(1-ratio);
m10 = 0;
m11 = ratio;
m12 = y0*(1-ratio);

% create transformation
trans = [m00 m01 m02; m10 m11 m12; 0 0 1];
