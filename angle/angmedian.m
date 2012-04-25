function med = angmedian(ang)
% ANGMEDIAN - Angular median
%
%   meg = angmedian(ang)
%
% Algorithm from Fisher, I think.

% Mercurial revision hash: $Revision$ $Date$
% See http://rcn.ccs.tulane.edu/index.php5/Tytell_Matlab
% Copyright (c) 2012, Eric Tytell <tytell at jhu dot edu>

ang = mod(ang,2*pi);
med = nanmedian(ang);

% average spacing between data points
space = nanmean(diff(sort(ang)));

% if we subtract the median, then take the median again, we should
% get zero.  If we don't, it means something has happened with the
% centering around zero.  The shift by pi is important in the mod,
% because otherwise we'd end up with things near 0 and near 2pi
ang2 = mod(ang-repmat(med,[size(ang,1) 1])+pi,2*pi) - pi;
dmed = nanmedian(ang2);

off = find(abs(dmed) > space);

% continue subtracting a median, recalculating the median and subtracting
% until we get something that is stable
reps = 0;
while (~isempty(off)),
    ang2(:,off) = mod(ang(:,off) - ...
                      repmat(med(:,off),[size(ang,1) 1])+pi,2*pi) - pi;
    dmed(off) = nanmedian(ang2(:,off));
    med(off) = med(off) + dmed(off);

    off = find(abs(dmed) > space);

    reps = reps+1;
    if (reps > 100),
        error('Too many iterations.');
    end;
end;


