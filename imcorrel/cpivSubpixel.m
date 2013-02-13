function frac = cpivSubpixel(val)

% remove negative values
m = min(val);
k = find(m <= 0);
if (~isempty(k)),
	m(k) = m(k) - 0.0001;		% everything should be slightly positive
	val(:,k) = val(:,k) - repmat(m(k), [3 1]);
end;

lval = log(val);

frac = (lval(1,:) - lval(3,:))./(2*(lval(1,:) + lval(3,:) - 2*lval(2,:)));
