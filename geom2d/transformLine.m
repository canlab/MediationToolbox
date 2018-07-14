function dest = transformLine(line, trans)
%TRANSFORMLINE tranform a line with an affine transform
%
%   LINE2 = transformLine(PT1, TRANS).
%   where LINE1 has the form [x0 y0 dx dy], and TRANS is a transformation
%   matrix, return the line transformed with affine transform TRANS. 
%
%   Format of TRANS can be one of :
%   [a b]   ,   [a b c] , or [a b c]
%   [d e]       [d e f]      [d e f]
%                            [0 0 1]
%
%   LINE2 = transformLine(LINES, TRANS) also wotk when LINES is a [N*4]
%   array of double. In this case, LINE2 has the same size as LINE.
%
%   See also :
%   transformPoint, translation, rotation, translation
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 06/04/2004.
%

% allocate memory
dest = zeros(size(line));

% compute position
dest(:,1) = line(:,1)*trans(1,1) + line(:,2)*trans(1,2);
dest(:,2) = line(:,1)*trans(2,1) + line(:,2)*trans(2,2);
dest(:,3) = line(:,3)*trans(1,1) + line(:,3)*trans(1,2);
dest(:,4) = line(:,4)*trans(2,1) + line(:,4)*trans(2,2);

% add translation vector, if exist
if size(trans, 2)>2
    dest(:,1) = dest(:,1)+trans(1,3);
    dest(:,2) = dest(:,2)+trans(2,3);
end