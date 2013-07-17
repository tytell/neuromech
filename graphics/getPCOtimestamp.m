function [dv,imnum] = getPCOtimestamp(I)

v = I(1,1:14);

if (any(v > 255))
    error('getpcotimestamp:notimestamp','No timestamp data');
end
v = uint8(v);

v1 = bitand(v,15);
v2 = bitshift(v,-4);
v = v2*10 + v1;

v = double(v);

imnum = v(4) + v(3)*100 + v(2)*100^2 + v(1)*100^3; 

yr = v(5)*100 + v(6);

month = v(7);
day = v(8);

h = v(9);
m = v(10);
s = v(11);

us1 = double(v(12:14));
us = us1(1)*10000 + us1(2)*100 + us1(3);

sec = double(s) + us/1e6;

dv = [double([yr month day h m]) sec];

