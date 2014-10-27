function [D] = show2DVec( Frame )
%Display located vector field 
D = create2DVec( Frame );
quiver( D.X, D.Y, D.U, D.V);
