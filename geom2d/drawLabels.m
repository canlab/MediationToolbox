function varargout = drawLabels(varargin)
%DRAWLABELS draw labels at specified positions
%   
%   DRAWLABELS(X, Y, LBL) draw labels LBL at position X and Y.
%   LBL can be either a string array, or a number array. In this case,
%   string are created by using sprintf function, with '%.2f' mask.
%
%   DRAWLABELS(POS, LBL) draw labels LBL at position specified by POS,
%   where POS is a N*2 int array.
%
%   DRAWLABELS(..., FORMAT, NUMBERS) create labels using sprintf function,
%   with the mask given by FORMAT (e. g. '%03d' or '5.3f'), and the
%   corresponding values.
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 15/12/2003.
%

if length(varargin)==1
    edge = varargin{1};
elseif length(varargin)==2
    p1 = varargin{1};
    p2 = varargin{2};
    edge = [p1(1) p1(2) p2(1)-p1(1) p2(2)-p1(2)];
end
    
h=line([edge(1) edge(1)+edge(3)], [edge(2) edge(2)+edge(4)]);
    

if nargout>0
    varargout{1}=h;
end