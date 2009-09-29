function [ind,sgn] = findzeros(y,varargin)
% function [ind,sgn] = findzeros(y, options...)
% Returns the index and sign (1 for ascending, -1 for descending) of zeros in y.

opt.numneighbors = 1;
opt.isstrict = false;
% opt.minzerodistance = 0;

opt = parsevarargin(opt,varargin,2);

if ((ndims(y) ~= 2) || all(size(y) ~= 1)),
    error('Can only operate on vectors');
end;

if (size(y,1) == 1)
    y = y';
    istranspose = true;
else
    istranspose = false;
end;

isasc = false(size(y));
isdesc = false(size(y));
isasc(1:end-1) = (y(1:end-1) < 0) & (y(2:end) >= 0);
isdesc(1:end-1) = (y(1:end-1) > 0) & (y(2:end) <= 0);

k = opt.numneighbors+1:length(y)-opt.numneighbors-1;
for off = 1:opt.numneighbors,
    isasc(k) = isasc(k) & (y(k-off) < 0) & (y(k+off+1) > 0);
    isdesc(k) = isdesc(k) & (y(k-off) > 0) & (y(k+off+1) < 0);
    
    if (opt.isstrict),
        isasc(k) = isasc(k) & (y(k-off) < y(k-off+1)) & (y(k+off+1) > y(k+off));
        isdesc(k) = isdesc(k) & (y(k-off) > y(k-off+1)) & (y(k+off+1) < y(k+off));
    end;
end;

iszero = isasc | isdesc;
ind = find(iszero);
sgn = zeros(size(y));
sgn(isasc) = 1;
sgn(isdesc) = -1;
sgn = sgn(iszero);

if (istranspose),
    ind = ind';
    sgn = sgn';
end;

    

    