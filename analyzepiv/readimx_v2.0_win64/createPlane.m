function [D]=createPlane( plane, scales )

%Create image data
nx = size(plane,1); 
ny = size(plane,2); 
D.X = double(1:nx)  * scales.X.Slope + scales.X.Offset ; 
D.Y = double(1:ny)  * scales.Y.Slope + scales.Y.Offset ;
D.I = double(plane) * scales.I.Slope + scales.I.Offset ;
