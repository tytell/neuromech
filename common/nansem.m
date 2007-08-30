function sem = nansem(x,dim)

if (nargin == 1)
    dim = [];
end;

sem = nanstd(x,[],dim) ./ sum(isfinite(x),dim);
