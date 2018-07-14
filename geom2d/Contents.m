% Geometry Toolbox
% Version 0.3 - 2005, November 7.
%
% 	Library to handle and visualize geometric primitives such as points,
% 	lines, circles and ellipses, polygons...
%
%   The goal is to provide a low-level library for manipulating geometrical
%   primitives, making easier the development of more complex geometric
%   algorithms. 
%
%
% Points :
% --------
%   angle3Points         - return oriented angle made by 3 points
%   centroid             - compute centroid (center of mass) of set of points
%   intersectEdges       - return intersections points of N edges in 2D
%   intersectLineEdge    - return intersection between a line and an edge
%   intersectLinePolygon - get intersection pts between a line and a polygon
%   intersectLines       - return intersection point of N lines in 2D
%   transformPoint       - tranform a point with an affine transform
%   pointOnLine          - create a point on a line at a given distance from line origin
%   polarPoint           - create a point from polar coordinates (rho + theta)
%   projPointOnLine      - return the projection of a point on a line
%   angleSort            - sort points of plane according to their angles
%
% Vectors :
% ---------
%   vecnorm              - compute norm of vector or of set of vectors
%   normalize            - normalize a vector
%
% Lines, edges and rays :
% -----------------------
%   createEdge           - create an edge between two points, or from a line
%   createLine           - create a line with various inputs.
%   createMedian         - create a median line
%   medianLine           - create a median line between two points
%   cartesianLine        - create a line with cartesian coefficients
%   bisector             - return the bisector of two lines, or 3 points
%   clipEdge             - clip an edge with a rectangular window
%   clipLineRect         - clip a line with a polygon
%   invertLine           - return same line but with opposite orientation
%   transformEdge        - tranform an edge with an affine transform
%   transformLine        - tranform a line with an affine transform
%   lineFit              - least mean square line regression
%   orthogonalLine       - create a line orthogonal to another one.
%   parallelLine         - create a line parallel to another one.
%
% Polygons :
% ----------
%   clipPolygon          - clip a polygon with a rectangular window
%   clipPolygonHP        - clip a polygon with Half-plane defined with directed line
%   convexification      - compute convexification of a polygon
%   medialAxisConvex     - compute medial axis of a convex polygon
%   polygonCentroid      - compute centroid (center of mass) of a polygon
%   polygonClipHP        - clip a polygon with Half-plane defined with directed line
%   polygonExpand        - 'expand' a polygon with a given distance
%   polygonArea          - compute area of a polygon
%   polygonLength        - compute perimeter of a polygon
%   readPolygon          - read a polygon stored in a file
%   steinerPoint         - compute steiner point (weighted centroid) of a polygon
%   steinerPolygon       - create a Steiner polygon from a set of vectors
%   supportFunction      - compute support function of a polygon
%
% Circles and Ellipses:
% ----------------------
%   createCircle         - create a circle from points
%   createDirectedCircle - create a directed circle
%   circleAsPolygon      - convert a circle into a series of points
%   ellipseAsPolygon     - convert an ellipse into a series of points
%   circleArcAsCurve     - convert a circle arc into a series of points
%   enclosingCircle      - find the minimum circle enclosing a set of points.
%
% Measurements:
% --------------------
%   distancePointEdge    - compute distance between a point and an edge
%   distancePointLine    - compute distance between a point and a line
%   distancePoints       - compute distance between two points
%   edgeAngle            - return angle of edge
%   edgeLength           - return length of an edge
%   lineAngle            - return angle between lines
%   linePosition         - return position of a point on a line
%   minDistance          - compute minimum distance between a point and a set of points
%   minDistancePoints    - compute minimal distance between several points
%   polygonNormalAngle   - compute normal angle at a vertex of the polygon
%
% Tests on shapes:
% --------------------
%   onEdge               - test if a point belongs to an edge
%   onLine               - test if a point belongs to a line
%   onRay                - test if a point belongs to a ray
%   onCircle             - test if a point is located on a given circle.
%   inCircle             - test if a point is located inside a given circle.
%   isLeftOriented       - check if a point is on the left side of a line
%
% Differential Geometry :
% -----------------------
%   polyfit2             - polynomial approximation of a curve
%   curveLength          - return length of a curve (a list of points)
%   curveCentroid        - compute centroid of a curve defined by a series of points
%   parametrize          - return a parametrization of a curve
%   curvature            - estimate curvature of a curve defined by points
%   surfaceCurvature     - compute curvature on a surface in a given direction 
%
% Other shapes :
% --------------
%   rectAsPolygon        - convert a (centered) rectangle into a series of points
%   crackPattern         - create a (bounded) crack pattern tessellation
%   crackPattern2        - create a (bounded) crack pattern tessellation
%   squareGrid           - generate equally spaces points in plane.
%   hexagonalGrid        - generate hexagonal grid of points in the plane.
%   triangleGrid         - generate triangular grid of points in the plane.
%
% Geometric transforms :
% ----------------------
%   cart2geod            - convert cartesian coordinates to geodesic coord.
%   geod2cart            - convert geodesic coordinates to cartesian coord.
%   homothecy            - create a homothecy as an affine transform
%   lineSymmetry         - create line symmetry as 2D affine transform
%   rotation             - return 3*3 matrix of a rotation
%   translation          - return 3*3 matrix of a translation
%   scaling              - return 3*3 matrix of a scale in 2 dimensions
%
%
% Drawing functions :
% -------------------
%   drawArrow            - draw an arrow on the current axis
%   drawCenteredEdge     - Draw an edge centered on a point
%   drawCircle           - draw a circle on the current axis
%   drawCircleArc        - draw a circle arc on the current axis
%   drawCurve            - Draw a curve specified by a list of points
%   drawEdge             - draw the edge given by 2 points
%   drawEllipse          - draw an ellipse on the current axis
%   drawEllipseArc       - draw an ellipse on the current axis
%   drawParabola         - draw a parabola on the current axis
%   drawLabels           - draw labels at specified positions
%   drawLine             - draw the line on the current axis
%   drawPoint            - draw the point on the axis.
%   drawPolygon          - draw a polygon specified by a list of points
%   drawRay              - draw a ray on the current axis
%   drawRect             - draw rectangle on the current axis
%   drawRect2            - draw centered rectangle on the current axis
%   drawShape            - draw various types of shapes (circles, polygons ...)
%   fillPolygon          - fill a polygon specified by a list of points
%
%
%   Credits:
%   * function 'enclosingCircle' rewritten from a file from Yazan Ahed
%       (yash78@gmail.com), available on Matlab File Exchange
%
%   -----
%
%   author : David Legland
%   INRA URPOI (Nantes) & MIA (Jouy-en-Josas)
%   david.legland@nantes.inra.fr
%   created the 07/11/2005.
%   Licensed under the terms of the LGPL, see the file "license.txt'

