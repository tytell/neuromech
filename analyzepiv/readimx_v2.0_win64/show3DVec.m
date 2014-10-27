function [D] = show3DVec( Frame )
%Display located vector field 
D = create3DVec( Frame );
quiver3(D.X,D.Y,D.Z,D.U,D.V,D.W,0);

