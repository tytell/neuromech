function [D]=showPlane( plane, scales )

D = createPlane( plane, scales );

%Display image data
imagesc(D.X,D.Y,D.I');
%Setup display
xlabel([scales.X.Description ' ' scales.X.Unit]);
ylabel([scales.Y.Description ' ' scales.Y.Unit]);
title ([scales.I.Description ' ' scales.I.Unit]);
if scales.Y.Slope >=0.0,
    axis ij;
else
    axis xy;
end
