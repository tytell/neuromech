function [D]=showVolume( planeCells, scales )

nz = size(planeCells,1); 
if nz < 1, error('invalid number of planes'); end;
nx = size(planeCells{1},1); 
ny = size(planeCells{1},2); 
%Create location data from scales
rx = double(1:nx)  * scales.X.Slope + scales.X.Offset ; 
ry = double(1:ny)  * scales.Y.Slope + scales.Y.Offset ;
rz = double(1:nz)  * scales.Z.Slope + scales.Z.Offset ;
% Display volume slizes
D = createVolume( planeCells, scales );
slice(D.X,D.Y,D.Z,D.I,[rx(nx/2)],[ry(ny/2)],rz(1:nz/4:nz)  );
