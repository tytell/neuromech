function y = medsmooth(x, n)

n2 = n/2;

wnd = repmat(1:n,[size(x,1) 1]) + repmat((0:size(x,1)-1)'-round(n2), [1 n]);
p = find((wnd >= 1) & (wnd <= size(x,1)));

val = zeros(size(wnd));
val(p) = x(wnd(p));

y = median(val,2);
