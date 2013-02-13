function DU = velderiv(x,y,u,v,varargin)
% function DU = velderiv(x,y,u,v,[w],[method],[sigma])
% takes derivatives of velocity vector (possibly stereo)

if ((nargin > 4) & isnumeric(varargin{5}) & ...
    (ndims(varargin{3}) == ndims(varargin{5})) & ...
    all(size(varargin{3}) == size(varargin{5}))),
    w = varargin{5};
    opts = varargin(6:end);
else
    w = [];
    opts = varargin(5:end);
end;

if isempty(w), 
    bWderiv = false; 
else, 
    bWderiv = true; 
end

if (length(opts) > 0),
    method = opts{1};
    if (length(opts) > 1),
        sigma = opts{2};
    end;
else
    method = 'central';
end;
method = lower(method);

sz = size(u);
if (ndims(u) > 2),
    szAll = sz;
    sz = sz(1:2);
end;
N = size(u,3);

if ((size(u,1) <= 2) | (size(u,2) <= 2)),
  error('Velocity matrix must be at least 3x3.');
end;

for i = 1:N,
    DU(i).dudx = repmat(NaN,sz);
    DU(i).dudy = repmat(NaN,sz);
    DU(i).dvdx = repmat(NaN,sz);
    DU(i).dvdy = repmat(NaN,sz);
    if (bWderiv),
        DU(i).dwdx = repmat(NaN,sz);
        DU(i).dwdy = repmat(NaN,sz);
    else
        DU(i).dwdx = [];
        DU(i).dwdy = [];
    end;

    if (size(x,3) > 1),
        x1 = x(:,:,i);
        y1 = y(:,:,i);
    else
        x1 = x;
        y1 = y;
    end;

    switch method,
     case 'central',
      DU(i).dudx(:,2:end-1) = (u(:,3:end,i) - u(:,1:end-2,i)) ./ ...
          (x1(:,3:end) - x1(:,1:end-2));
      DU(i).dudx(:,[1 end]) = (u(:,[2 end],i) - u(:,[1 end-1],i)) ./ ...
          (x1(:,[2 end]) - x1(:,[1 end-1]));
      DU(i).dudy(2:end-1,:) = (u(3:end,:,i) - u(1:end-2,:,i)) ./ ...
          (y1(3:end,:) - y1(1:end-2,:));
      DU(i).dudy([1 end],:) = (u([2 end],:,i) - u([1 end-1],:,i)) ./ ...
          (y1([2 end],:) - y1([1 end-1],:));

      DU(i).dvdx(:,2:end-1) = (v(:,3:end,i) - v(:,1:end-2,i)) ./ ...
          (x1(:,3:end) - x1(:,1:end-2));
      DU(i).dvdx(:,[1 end]) = (v(:,[2 end],i) - v(:,[1 end-1],i)) ./ ...
          (x1(:,[2 end]) - x1(:,[1 end-1]));
      DU(i).dvdy(2:end-1,:) = (v(3:end,:,i) - v(1:end-2,:,i)) ./ ...
          (y1(3:end,:) - y1(1:end-2,:));
      DU(i).dvdy([1 end],:) = (v([2 end],:,i) - v([1 end-1],:,i)) ./ ...
          (y1([2 end],:) - y1([1 end-1],:));
      if bWderiv
          DU(i).dwdx(:,2:end-1) = (w(:,3:end,i) - w(:,1:end-2,i)) ./ ...
              (x1(:,3:end) - x1(:,1:end-2));
          DU(i).dwdx(:,[1 end]) = (w(:,[2 end],i) - w(:,[1 end-1],i)) ./ ...
              (x1(:,[2 end]) - x1(:,[1 end-1]));
          DU(i).dwdy(2:end-1,:) = (w(3:end,:,i) - w(1:end-2,:,i)) ./ ...
              (y1(3:end,:) - y1(1:end-2,:));
          DU(i).dwdy([1 end],:) = (w([2 end],:,i) - w([1 end-1],:,i)) ./ ...
              (y1([2 end],:) - y1([1 end-1],:));
      end;

     case 'nogueira',
      DU(i).dudx = nogueiraderiv(x,u(:,:,i), sigma);
      DU(i).dudy = nogueiraderiv(y',u(:,:,i)', sigma)';
      DU(i).dvdx = nogueiraderiv(x,v(:,:,i), sigma);
      DU(i).dvdy = nogueiraderiv(y',v(:,:,i)', sigma)';
      if bWderiv
          DU(i).dwdx = nogueiraderiv(x,w(:,:,i), sigma);
          DU(i).dwdy = nogueiraderiv(y',w(:,:,i)', sigma)';
      end
    end;  % case
end;


function d = nogueiraderiv(x,y, sigma)

Fj = 1/1900*[-9 -57 404 -1551 0 1551 -404 57 9]';
Fk = 1/1900*[124 -93 -398 -371 0 371 398 93 -124]';
Fsign = [-4 0 1 14 0 -14 -1 0 4; ...
         4.5 7.9 -9 23.7 0 -23.7 9 7.9 -4.5; ...
         1 0 -8 12 0 -12 8 0 -1; ...
         -1 1 -2 5 0 -5 2 -1 1]';

A = im2col(y,[1 9],'sliding');
d = sum(A.*repmat(Fj,[1 size(A,2)]));
lo = find(d < 3.5*sigma);
hi = find(d >= 3.5*sigma);
sgn(1,:) = sum(A(:,hi).*repmat(Fsign(:,1),[1 length(hi)]));
sgn(2,:) = sum(A(:,hi).*repmat(Fsign(:,2),[1 length(hi)]));
sgn(3,:) = sum(A(:,hi).*repmat(Fsign(:,3),[1 length(hi)]));
sgn(4,:) = sum(A(:,hi).*repmat(Fsign(:,4),[1 length(hi)]));
sbad = hi(find(all(sign(sgn) == repmat(sign(sgn(1,:)),[4 1]))));		

d([lo sbad]) = sum(A(:,[lo sbad]).*repmat(Fk,[1 length(lo)+length(sbad)]));
d = d/(x(1,2)-x(1));

d = col2im(d,[1 9],size(y),'sliding');
d = [repmat(NaN,[size(y,1) 4]) d repmat(NaN,[size(y,1) 4])];
