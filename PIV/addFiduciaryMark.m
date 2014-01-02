function addFiduciaryMark(imname)

if (nargin == 0) || isempty(imname)
    [fn,pn] = uigetfile('*.tif','Select calibration image');
    imname = fullfile(pn,fn);
end

[pn,fn,ext] = fileparts(imname);

I = imread(imname);
if (isa(I,'uint16')) && (max(I(:)) <= 2^12)
    % probably a 12 bit image
    I = im2uint8(I*16);
end

figure(1);
imagesc(I);
colormap(gray);

zoom on;
input('Zoom to center point and hit return');
zoom off;

disp('Click center point');
[j,i] = ginputb(1);

done = false;
while (~done)
    sz = input('Size? ');
    ln = round(max(sz)/10);
    if ln < 1
        ln = 1;
    end

    I2 = I;
    i = round(i);
    j = round(j);
    val = I(i,j);

    if (numel(sz) == 1)
        sz = [sz sz];
    end
    sz2 = round(sz/2);
    I2(i-ln:i+ln, j-sz2(2):j+sz2(2)) = val;
    I2(i-sz2(1):i+sz2(1), j-ln:j+ln) = val;

    figure(2);
    imagesc(I2);
    
    done = inputyn('Done? ');
end


[fn,pn] = uiputfile('*.tif','Choose output file',[fn '-mod' ext]);
imwrite(I2, fullfile(pn,fn));



