function [D] = create2DVec( Frame )

Components = Frame.Components;
frameInfo = MakeFrameInfo(Frame);
%Compile vector data
if (frameInfo.hasChoices),
	D.C = Components{ frameInfo.best }.Planes{1};
	D.U = zeros(size(D.C));
	D.V = zeros(size(D.C));
	for i=0:4,
		mask = (D.C==i);
        j = min(i,3);
        choiceU = frameInfo.choices(j*3+1);
        choiceV = frameInfo.choices(j*3+2);
		D.U(mask) = Components{ choiceU }.Planes{1}(mask);
		D.V(mask) = Components{ choiceV }.Planes{1}(mask);
	end
else
	D.U = Components{frameInfo.choices(1)}.Planes{1}; 
	D.V = Components{frameInfo.choices(2)}.Planes{1}; 
end
Grids  = frameInfo.Grids;
scaleX = frameInfo.Scales.X; 
scaleY = frameInfo.Scales.Y; 
scaleI = frameInfo.Scales.I;
D.U = double(D.U)*scaleI.Slope + scaleI.Offset;
D.V = double(D.V)*scaleI.Slope + scaleI.Offset;
if scaleY.Slope < 0.0, D.V = -D.V; end
%Compile location data
Rx = (1:size(D.U,1))-0.5; 
Ry = (1:size(D.U,2))-0.5;
[D.X,D.Y] = ndgrid(Rx*Grids.X*scaleX.Slope + scaleX.Offset, Ry*Grids.Y*scaleY.Slope + scaleY.Offset);
