function line = invertLine(var)
%INVERTLINE return same line but with opposite orientation
%
%   iLine = invertLine(line) return the opposite line.
%   'line' has the format [x0 y0 dx dy], then 'iLine' will have following
%   parameters : [x0 y0 -dx -dy].
%
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 20/01/2004.
%

line = 0;    

if size(var, 1)==1
    % only one line in a single array
    line = [var(1) var(2) -var(3) -var(4)];
else
    % several lines in a single array
    n = size(var, 1);
    line(1:n, 1) = var(1:n, 1);
    line(1:n, 2) = var(1:n, 2);
    line(1:n, 3) = -var(1:n, 3);
    line(1:n, 4) = -var(1:n, 4);
end
