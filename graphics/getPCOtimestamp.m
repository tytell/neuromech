function [dv,imnum,raw] = getPCOtimestamp(I)

v = I(1,1:14);

if (any(v > 255))
    error('No timestamp data');
end
v = uint8(v);

imnum1 = v(1:4);
imnum = typecast(imnum1(end:-1:1),'uint32');

yr1 = v(5:6);
yr = typecast(yr1([2 1]),'uint16');

month = v(7);
day = v(8);

h = v(9);
m = v(10);
s = v(11);

us1 = double(v(12:14));
us = us1(1)*10000 + us1(2)*100 + us1(3);

sec = double(s) + us/1e6;

dv = [double([yr month day h m]) sec];
raw = v;
