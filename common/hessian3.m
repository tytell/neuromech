function H = hessian3(x,y,z, V, method)
% function H = hessian3(x,y,z, V, method)
% Calculates a 3D matrix of second derivatives using a 5th order
% central difference algorithm.
% Copyright (c) 2004 Eric Tytell

if (nargin == 4),
  method = 'difference';
end;

switch lower(method),
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

  [j,i,k] = meshgrid(1:size(V,2), 1:size(V,1), 1:size(V,3));

  off = repmat(inda,[1 size(V)]);
  aa1 = repmat(a1,[1 size(V)]);
  aa2 = repmat(a2,[1 size(V)]);
  ii = repmat(shiftdim(i,-1), [2*ord+1 1 1 1]);
  jj = repmat(shiftdim(j,-1), [2*ord+1 1 1 1]);
  kk = repmat(shiftdim(k,-1), [2*ord+1 1 1 1]);

  H = repmat(NaN,[3 3 size(V)]);

  m1 = ord+1:size(V,1)-ord;
  ind = sub2ind(size(V), ii(:,m1,:,:)+off(:,m1,:,:),...
                jj(:,m1,:,:), kk(:,m1,:,:));
  H(2,2,m1,:,:) = 1/dy^2 * squeeze(sum(V(ind).*aa2(:,m1,:,:)));

  H(1,2,m1,:,:) = 1/dy * squeeze(sum(V(ind).*aa1(:,m1,:,:)));

  m2 = ord+1:size(V,2)-ord;
  ind = sub2ind(size(V), ii(:,:,m2,:),...
                jj(:,:,m2,:)+off(:,:,m2,:), kk(:,:,m2,:));
  H(1,1,:,m2,:) = 1/dx^2 * squeeze(sum(V(ind).*aa2(:,:,m2,:)));

  H(1,2,:,m2,:) = 1/dx * squeeze(sum( ...
      reshape(H(1,2,ind),size(ind)).*aa1(:,:,m2,:)));
  H(1,2,:,[1:ord end-ord+1:end],:) = NaN;
  H(2,1,:,:,:) = H(1,2,:,:,:);

  if (size(V,3) >= 2*ord+1),
    H(1,3,:,m2,:) = 1/dx * squeeze(sum(V(ind).*aa1(:,:,m2,:)));

    m3 = ord+1:size(V,3)-ord;
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
  if (length(z) > 1),
    sp = spapi({6,6,6}, {y,x,z}, V);
    H(1,1,:,:,:) = fnval(fnder(sp,[0 2 0]),{y,x,z});
    H(1,2,:,:,:) = fnval(fnder(sp,[1 1 0]),{y,x,z});
    H(2,1,:,:,:) = H(1,2,:,:,:);
    H(1,3,:,:,:) = fnval(fnder(sp,[0 1 1]),{y,x,z});
    H(3,1,:,:,:) = H(1,3,:,:,:);
    H(2,2,:,:,:) = fnval(fnder(sp,[2 0 0]),{y,x,z});
    H(2,3,:,:,:) = fnval(fnder(sp,[1 0 1]),{y,x,z});
    H(3,2,:,:,:) = H(2,3,:,:,:);
    H(3,3,:,:,:) = fnval(fnder(sp,[0 0 2]),{y,x,z});
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




