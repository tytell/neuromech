function [x,y,I2,ix,iy] = digitizeGraph(I)

if (nargin == 0),
    [fn,pn] = uigetfile({'*.jpg;*.tif;*.png;*.bmp','Image files'}, ...
        'Graph image file');
    fname = fullfile(pn,fn);
    
    I = imread(fname);
end;

imshow6(I,'n');

box = msgbox('Click origin','Digitize graph','help','non-modal');
[ox,oy] = ginput(1);
if (ishandle(box)),
    delete(box);
end;

box = msgbox('Click x axis limit','Digitize graph','help','non-modal');
[xaxx,xaxy] = ginput(1);
if (ishandle(box)),
    delete(box);
end;

box = msgbox('Click y axis limit','Digitize graph','help','non-modal');
[yaxx,yaxy] = ginput(1);
if (ishandle(box)),
    delete(box);
end;

lims = inputdlg({'X axis minimum:', 'X axis maximum:',...
    'Y axis minimum:','Y axis maximum:'},'Axes limits',1, ...
    {'0','','0',''});
lims = cellfun(@str2double,lims);

tform = cp2tform([ox oy; xaxx xaxy; yaxx yaxy],...
    [lims(1) lims(3); lims(2) lims(3); lims(1) lims(4)], 'affine');

[I2,ix,iy] = imtransform(I, tform, 'Size',[size(I,1) size(I,2)]);

imshow6(ix,iy,I2,'n');
axis normal xy tight;

box = msgbox('Click points to digitize','Digitize graph','help','non-modal');
[x,y] = ginput;
if (ishandle(box)),
    delete(box);
end;
