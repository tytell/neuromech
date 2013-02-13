function mad = angmeanavgdev(ang)

ang = mod(ang,2*pi);
med = angmedian(ang);

mad = pi - nanmean(abs(pi - abs(ang - repmat(med,[size(ang,1) 1]))));
