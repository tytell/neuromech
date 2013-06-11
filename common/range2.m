function r = range2(x,dim)

if nargin == 1
    if (isrow(x))
        dim = 2;
    else
        dim = 1;
    end
end

ord = [dim 1:dim-1 dim+1:ndims(x)];
x = permute(x,ord);

sz = size(x);
r = NaN([1 sz(2:end)]);

good = isfinite(x);

n = prod(sz(2:end));
for i = 1:n
    if (any(good(:,i)))
        a = min(x(good(:,i),i));
        b = max(x(good(:,i),i));
        r(1,i) = b - a;
    end
end

r = ipermute(r,ord);

    
    

