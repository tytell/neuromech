function [frameInfo]=MakeFrameInfo(Frame)

Names = Frame.ComponentNames;
frameInfo.choices = zeros(12,1);
for i=1:size(Names,1),
	if strcmp(Names{i},'U0')==1, frameInfo.choices( 1) = i; end; 
	if strcmp(Names{i},'V0')==1, frameInfo.choices( 2) = i; end; 
	if strcmp(Names{i},'W0')==1, frameInfo.choices( 3) = i; end; 
	if strcmp(Names{i},'U1')==1, frameInfo.choices( 4) = i; end; 
	if strcmp(Names{i},'V1')==1, frameInfo.choices( 5) = i; end; 
	if strcmp(Names{i},'W1')==1, frameInfo.choices( 6) = i; end; 
	if strcmp(Names{i},'U2')==1, frameInfo.choices( 7) = i; end; 
	if strcmp(Names{i},'V2')==1, frameInfo.choices( 8) = i; end; 
	if strcmp(Names{i},'W2')==1, frameInfo.choices( 9) = i; end; 
	if strcmp(Names{i},'U3')==1, frameInfo.choices(10) = i; end; 
	if strcmp(Names{i},'V3')==1, frameInfo.choices(11) = i; end; 
	if strcmp(Names{i},'W3')==1, frameInfo.choices(12) = i; end; 
	if strcmp(Names{i},'ACTIVE_CHOICE')== 1, frameInfo.best   = i; end;
	if strcmp(Names{i},'ENABLED')      == 1, frameInfo.enable = i; end; 
	if strcmp(Names{i},'MASK')         == 1, frameInfo.mask   = i; end; 
end
frameInfo.Grids  = Frame.Grids;
frameInfo.Scales.X = Frame.Scales.X;
frameInfo.Scales.Y = Frame.Scales.Y;
frameInfo.Scales.Z = Frame.Scales.Z;
frameInfo.Scales.I = Frame.Scales.I;
frameInfo.is3D = (frameInfo.choices(3) > 0);	
frameInfo.hasChoices = (frameInfo.choices(4) > 0);