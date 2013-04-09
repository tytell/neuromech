function m = nanmedian2(x,dim)

ord = [dim 1:dim-1 dim+1:ndims(x)];
x = permute(x,ord);

sz = size(x);
m = zeros([1 sz(2:end)]);

good = isfinite(x);

n = prod(sz(2:end));
for i = 1:n
    m(1,i) = median(x(good(:,i),i));
end

m = ipermute(m,ord);

    
    

