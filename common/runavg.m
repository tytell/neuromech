function M = runavg(m, n)
% function M = runavg(m, n)

n2 = round(n/2);
len = length(m);
M = zeros(size(m));

for i = 1:n2,
   a = mean(m(1:i+n2));
   M(i) = a;
end;

for i = n2+1:len-n2,
    a = mean(m(i-n2:i+n2));
    M(i) = a; 
end;
   
for i = len-n2:len,
    a = mean(m(i-n2:len));
    M(i) = a;
end;

   	