function H = hessian3(x,y,z, V)

if (length(z) > 1),
  sp = spapi({6,6,6}, {x,y,z}, V);
  H(:,:,:,1,1) = fnval(fnder(sp,[2 0 0]),{x,y,z});
  H(:,:,:,1,2) = fnval(fnder(sp,[1 1 0]),{x,y,z});
  H(:,:,:,2,1) = H(:,:,:,1,2);
  H(:,:,:,1,3) = fnval(fnder(sp,[1 0 1]),{x,y,z});
  H(:,:,:,3,1) = H(:,:,:,1,3);
  H(:,:,:,2,2) = fnval(fnder(sp,[0 2 0]),{x,y,z});
  H(:,:,:,2,3) = fnval(fnder(sp,[0 1 1]),{x,y,z});
  H(:,:,:,3,2) = H(:,:,:,2,3);
  H(:,:,:,3,3) = fnval(fnder(sp,[0 0 2]),{x,y,z});
else
  sp = spapi({6,6}, {x,y}, V);
  H(:,:,:,1,1) = fnval(fnder(sp,[2 0]),{x,y});
  H(:,:,:,1,2) = fnval(fnder(sp,[1 1]),{x,y});
  H(:,:,:,2,1) = H(:,:,:,1,2);
  H(:,:,:,1,3) = 0;
  H(:,:,:,3,1) = 0;
  H(:,:,:,2,2) = fnval(fnder(sp,[0 2]),{x,y});
  H(:,:,:,2,3) = 0;
  H(:,:,:,3,2) = 0;
  H(:,:,:,3,3) = 0;
end;







