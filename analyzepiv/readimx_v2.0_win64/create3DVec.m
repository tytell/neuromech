function [D] = create3DVec( Frame )

Components = Frame.Components;
frameInfo = MakeFrameInfo(Frame);
%Compile best vector data
if (frameInfo.hasChoices),
	D.C = Components{ frameInfo.best }.Planes{1};
	D.U = zeros(size(D.C));
	D.V = zeros(size(D.C));
	D.W = zeros(size(D.C));
	for i=0:5,
		mask = (D.C==i);
        j = min(i,3);
        choiceU = frameInfo.choices(j*3+1);
        choiceV = frameInfo.choices(j*3+2);
        choiceW = frameInfo.choices(j*3+3);
		D.U(mask) = Components{ choiceU }.Planes{1}(mask);
		D.V(mask) = Components{ choiceV }.Planes{1}(mask);
		D.W(mask) = Components{ choiceW }.Planes{1}(mask);
	end
else
	D.U = Components{frameInfo.choices(1)}.Planes{1}; 
	D.V = Components{frameInfo.choices(2)}.Planes{1}; 
	D.W = Components{frameInfo.choices(3)}.Planes{1}; 
end
scaleI = frameInfo.Scales.I;
D.U = double(D.U)*scaleI.Slope + scaleI.Offset;
D.V = double(D.V)*scaleI.Slope + scaleI.Offset;
D.W = double(D.W)*scaleI.Slope + scaleI.Offset;
scaleY = frameInfo.Scales.Y; 
if scaleY.Slope < 0.0, D.V = -D.V; end
%Complile location data
Grids  = frameInfo.Grids;
scaleX = frameInfo.Scales.X; scaleZ = frameInfo.Scales.Z;
Rx = (1:size(D.U,1))-0.5; 
Ry = (1:size(D.U,2))-0.5; 
Rz = (1:size(D.U,3))-1.0;
[D.X,D.Y,D.Z] = ndgrid( Rx*Grids.X*scaleX.Slope + scaleX.Offset,...
                        Ry*Grids.Y*scaleY.Slope + scaleY.Offset,... 
				        Rz        *scaleZ.Slope + scaleZ.Offset);