function varargout = hessian3(x,y,z, V, method, ind1,ind2,ind3)
% function H = hessian3(x,y,z, V, method)
% Calculates a 3D matrix of second derivatives using a 5th order
% central difference algorithm.
% Copyright (c) 2004 Eric Tytell

if (nargin == 4),
  method = 'difference';
end;

switch lower(method),
 case 'ndifference',
  if ((nargin <= 6) | ((nargin == 8) & isempty(ind2))),
      if (exist('ind1') & ~isempty(ind1)),
          N = floor(ind1/2);
          if (2*N+1 > min(size(V))),
              error('Order is too large');
          end;
      else
          N = floor(min(size(V)-1)/2);
      end;
  
      %NB: I have no idea what this section here is supposed to do, and the
      %mfactorial program seems to have disappeared...  (ET 7/28/06)
      k = (-N:N)';
      CNk = mfactorial(N)^2./(mfactorial(N-k).*mfactorial(N+k));
      ctr = N+1;
      n = length(k);

      nonzero = [1:ctr-1 ctr+1:length(k)];
      first(nonzero,1) = (-1).^(k(nonzero)+1) .* 1./k(nonzero) .* ...
          CNk(nonzero);
      first(ctr) = 0;
      second(nonzero,1) = (-1).^(k(nonzero)+1) .* 2./k(nonzero).^2 .* ...
          CNk(nonzero);
      second(ctr) = -2*sum(second(ctr+1:end));
  else
      first = ind2;
      second = ind3;
      N = floor(length(first)/2);
      k = (-N:N)';
  end;

  mid = ceil(size(V)/2);
  dx = x(2)-x(1);
  dy = y(2)-y(1);
  dz = z(2)-z(1);

  H(1,1) = 1/dx^2 * sum(second.*V(mid(1),mid(2)+k,mid(3))');
  H(2,2) = 1/dy^2 * sum(second.*V(mid(1)+k,mid(2),mid(3)));
  H(3,3) = 1/dz^2 * sum(second.*squeeze(V(mid(1),mid(2),mid(3)+k)));

  cross = first*first';
  H(1,2) = 1/(dx*dy) * sum(sum(cross.*V(mid(1)+k,mid(2)+k,mid(3))));
  H(1,3) = 1/(dx*dz) * sum(sum(permute(cross,[3 1 2]) .* ...
                               V(mid(1),mid(2)+k,mid(3)+k)));
  H(2,3) = 1/(dy*dz) * sum(sum(permute(cross,[1 3 2]) .* ...
                               V(mid(1)+k,mid(2),mid(3)+k)));
  H(2,1) = H(1,2);
  H(3,1) = H(1,3);
  H(3,2) = H(2,3);

  if (nargout == 3),
      varargout{2} = first;
      varargout{3} = second;
  end;
  
  
 case 'difference',
  a1 = [-1 9 -45 0 45 -9 1]'/60;
  a2 = [2 -27 270 -490 270 -27 2]'/180;

  ord = 3;
  inda = (-ord:ord)';

  dx = x(2)-x(1);
  dy = y(2)-y(1);
  if (length(z) > 1),
    dz = z(2)-z(1);
  else
    dz = 0;
  end;

  if (exist('ind1')),
      error('This option doesn''t work.');

      i = ind1(1)-ord:ind1(end)+ord;
      j = ind2(1)-ord:ind2(end)+ord;
      k = ind3(1)-ord:ind3(end)+ord;
      i = i((i >= 1) & (i <= size(V,1)));
      j = j((j >= 1) & (j <= size(V,2)));
      k = k((k >= 1) & (k <= size(V,3)));

      [j,i,k] = meshgrid(j,i,k);

      m1 = ind1-ind1(1)+ord+1;
      m2 = ind2-ind2(1)+ord+1;
      m3 = ind3-ind3(1)+ord+1;
  else
      [j,i,k] = meshgrid(1:size(V,2), 1:size(V,1), 1:size(V,3));
      m1 = ord+1:size(V,1)-ord;
      m2 = ord+1:size(V,2)-ord;
      m3 = ord+1:size(V,3)-ord;
  end;

  off = repmat(inda,[1 size(i)]);
  aa1 = repmat(a1,[1 size(i)]);
  aa2 = repmat(a2,[1 size(i)]);
  ii = repmat(shiftdim(i,-1), [2*ord+1 1 1 1]);
  jj = repmat(shiftdim(j,-1), [2*ord+1 1 1 1]);
  kk = repmat(shiftdim(k,-1), [2*ord+1 1 1 1]);

  H = repmat(NaN,[3 3 size(i)]);

  ind = sub2ind(size(V), ii(:,m1,:,:)+off(:,m1,:,:),...
                jj(:,m1,:,:), kk(:,m1,:,:));
  H(2,2,m1,:,:) = 1/dy^2 * squeeze(sum(V(ind).*aa2(:,m1,:,:)));

  H(1,2,m1,:,:) = 1/dy * squeeze(sum(V(ind).*aa1(:,m1,:,:)));

  ind = sub2ind(size(V), ii(:,:,m2,:),...
                jj(:,:,m2,:)+off(:,:,m2,:), kk(:,:,m2,:));
  H(1,1,:,m2,:) = 1/dx^2 * squeeze(sum(V(ind).*aa2(:,:,m2,:)));

  H(1,2,:,m2,:) = 1/dx * squeeze(sum( ...
      reshape(H(1,2,ind),size(ind)).*aa1(:,:,m2,:)));
  H(1,2,:,[1:ord end-ord+1:end],:) = NaN;
  H(2,1,:,:,:) = H(1,2,:,:,:);

  if (size(V,3) >= 2*ord+1),
    H(1,3,:,m2,:) = 1/dx * squeeze(sum(V(ind).*aa1(:,:,m2,:)));

    ind = sub2ind(size(V), ii(:,:,:,m3), ...
                  jj(:,:,:,m3), kk(:,:,:,m3)+off(:,:,:,m3));
    H(3,3,:,:,m3) = 1/dz^2 * squeeze(sum(V(ind).*aa2(:,:,:,m3)));

    H(1,3,:,:,m3) = 1/dz * squeeze(sum( ...
        reshape(H(1,3,ind),size(ind)).*aa1(:,:,:,m3)));
    H(1,3,:,:,[1:ord end-ord+1:end]) = NaN;
    H(3,1,:,:,:) = H(1,3,:,:,:);

    H(2,3,:,:,m3) = 1/dz * squeeze(sum(V(ind).*aa1(:,:,:,m3)));
    ind = sub2ind(size(V), ii(:,m1,:,:)+off(:,m1,:,:),...
                  jj(:,m1,:,:), kk(:,m1,:,:));
    H(2,3,m1,:,:) = 1/dy * squeeze(sum( ...
        reshape(H(2,3,ind),size(ind)).*aa1(:,m1,:,:)));
    H(2,3,[1:ord end-ord+1:end],:,:) = NaN;
    H(3,2,:,:,:) = H(2,3,:,:,:);
  else
    H = H(1:2,1:2,:,:);
  end;
  
 case 'spline',
  ord = 6;
  noff = ceil(ord/2)+2;

  if (exist('ind1')),
      i = ind1(1)-noff:ind1(end)+noff;
      j = ind2(1)-noff:ind2(end)+noff;
      k = ind3(1)-noff:ind3(end)+noff;
      i = i((i >= 1) & (i <= size(V,1)));
      j = j((j >= 1) & (j <= size(V,2)));
      k = k((k >= 1) & (k <= size(V,3)));
  else
      i = 1:size(V,1);
      j = 1:size(V,2);
      k = 1:size(V,3);
      ind1 = i;
      ind2 = j;
      ind3 = k;
  end;

  if (length(z) > 1),
    sp = spapi({ord,ord,ord}, {y(i),x(j),z(k)}, V(i,j,k));
    H(1,1,:,:,:) = fnval(fnder(sp,[0 2 0]),{y(ind1),x(ind2),z(ind3)});
    H(1,2,:,:,:) = fnval(fnder(sp,[1 1 0]),{y(ind1),x(ind2),z(ind3)});
    H(2,1,:,:,:) = H(1,2,:,:,:);
    H(1,3,:,:,:) = fnval(fnder(sp,[0 1 1]),{y(ind1),x(ind2),z(ind3)});
    H(3,1,:,:,:) = H(1,3,:,:,:);
    H(2,2,:,:,:) = fnval(fnder(sp,[2 0 0]),{y(ind1),x(ind2),z(ind3)});
    H(2,3,:,:,:) = fnval(fnder(sp,[1 0 1]),{y(ind1),x(ind2),z(ind3)});
    H(3,2,:,:,:) = H(2,3,:,:,:);
    H(3,3,:,:,:) = fnval(fnder(sp,[0 0 2]),{y(ind1),x(ind2),z(ind3)});
  else
    sp = spapi({6,6}, {x,y}, V);
    H(1,1,:,:,:) = fnval(fnder(sp,[0 2]),{y,x});
    H(1,2,:,:,:) = fnval(fnder(sp,[1 1]),{y,x});
    H(2,1,:,:,:) = H(1,2,:,:,:);
    H(1,3,:,:,:) = 0;
    H(3,1,:,:,:) = 0;
    H(2,2,:,:,:) = fnval(fnder(sp,[2 0]),{y,x});
    H(2,3,:,:,:) = 0;
    H(3,2,:,:,:) = 0;
    H(3,3,:,:,:) = 0;
  end;
end;

varargout{1} = H;




