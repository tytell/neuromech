function M = runavg(m, n, dim)
% function M = runavg(m, n, dim)

if (nargin == 2)
    dim = 1;
end;

if (size(m,dim) < n),
    error('Too few data points');
end;

%permute and reshape m so that it has the dimension of interest along the
%rows and everything else flatten along columns
pmt = [dim 1:dim-1 dim+1:ndims(m)];
m = permute(m,pmt);
szm = size(m);
m = reshape(m,[szm(1) prod(szm(2:end))]);

%get the sizes
n2 = floor(n/2);
n = n2*2 + 1;
len = size(m,1);

%set up the output matrix
M = zeros(size(m));

%get the indices of most of the values
ind = repmat((-n2:n2)',[1 len]) + repmat(1:len,[n 1]);
%deal with the beginning separately
ind1 = ind(:,1:n2);
good1 = ind1(:,1:n2) >= 1;
%and the end
ind2 = ind(:,len-n2+1:len);
good2 = ind(:,len-n2+1:len) <= len;

%the middle bit
k = n2+1:len-n2;
ind = ind(:,k);

%number of values
N1 = sum(good1);
N2 = sum(good2);

%run through the columns
mm1 = zeros(n,n2);
mm2 = zeros(n,n2);

for i = 1:size(m,2),
    mm = reshape(m(ind,i),size(ind));
    mm1(good1) = m(ind1(good1),i);
    mm2(good2) = m(ind2(good2),i);
    
    M(1:n2,i) = sum(mm1)./N1;
    M(k,i) = sum(mm)/n;
    M(len-n2+1:len,i) = sum(mm2)./N2;
end;

%and set M back to the same shape and size as m
M = reshape(M,szm);
M = ipermute(M,pmt);

   