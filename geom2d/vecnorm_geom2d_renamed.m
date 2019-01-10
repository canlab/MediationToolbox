function n = vecnorm(v, varargin)
%VECNORM compute norm of vector or of set of vectors
%
%   n = vecnorm(V);
%   return euclidean norm of vector V.
%
%   n = vecnorm(V, N);
%   specify the norm to use. N can be any value greater than 0. 
%   N=1 -> city lock norm
%   N=2 -> euclidean norm
%   N=inf -> compute max coord.
%
%   When V is a MxN array, compute norm for each vector of the array.
%   Vector are given as rows. Result is then a [M*1] array.
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 21/02/2005.
%

%   HISTORY
%   02/05/2006  manage several norms

dim = size(v);

d = 2;
if length(varargin)>0
    d = varargin{1};
end

if d==1
    if dim(1)==1 || dim(2)==1
        n = sum(abs(v));
    else
        n = sum(abs(v), 2);
    end
elseif d==2
    if dim(1)==1 || dim(2)==1
        n = sqrt(sum(v.*v));
    else
        n = sqrt(sum(v.*v, 2));
    end
elseif d==inf
    if dim(1)==1 || dim(2)==1
        n = max(v);
    else
        n = max(v, 2);
    end
else
    if dim(1)==1 || dim(2)==1
        n = power(sum(power(v, d)), 1/d);
    else
        n = power(sum(power(v, d), 2), 1/d);
    end
end    