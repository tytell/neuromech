function [D]=createVolume( planeCells, scales )

nz = size(planeCells,1); 
if nz < 1, error('invalid number of planes'); end;
nx = size(planeCells{1},1); 
ny = size(planeCells{1},2); 
%Create location data from scales
rx = double(1:nx) * scales.X.Slope + scales.X.Offset ; 
ry = double(1:ny) * scales.Y.Slope + scales.Y.Offset ;
rz = double(1:nz) * scales.Z.Slope + scales.Z.Offset ;
[D.X,D.Y,D.Z] = meshgrid(rx,ry,rz);
%Create volume data from planes
D.I = zeros(nx,ny,nz);
for i=1:nz, 
    D.I(:,:,i) = planeCells{i}* scales.I.Slope + scales.I.Offset; 
end;
