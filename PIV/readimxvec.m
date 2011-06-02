function varargout = readimxvec(file,choice)
% function [x,y,u,v,w?,units] = readimxvec(file)

if (nargin == 1),
    choice = [];
end;

A = readimx(file);

nx = A.Nx;
ny = A.Ny;

rx = 1:nx;
ry = 1:ny;

%set up the coordinate system
x = (rx-1)*A.Grid*A.ScaleX(1) + A.ScaleX(2);
y = (ry-1)*A.Grid*A.ScaleY(1) + A.ScaleY(2);

switch A.IType,
 case 0,
  error('Cannot read image files.');

 case {1,3,5},
  if (A.IType == 5),
      nvec = 3;                         % 3 component vectors
  else
      nvec = 2;                         % 2 component vectors
  end;

  npages = size(A.Data,2)/ny;
  %error if we can't divide A.Data into a discrete number of pages
  if (mod(npages,1) ~= 0),
      error('Non-integer number of pages in Data structure.');
  end;

  data = double(reshape(A.Data,[nx ny npages]));

  %first page in data is an 1-based index into successive pages for which
  %vector DaVis chose.  Last pages are postprocessed vectors (given indices
  %4 or 5).  But since there are a max of 4*nvec pages, we have to
  %truncate at pageind == 3.  Don't ask me why
  if (~isempty(choice)),
      pageind = repmat(choice,[size(data,1) size(data,2)]);
  else
      pageind = data(:,:,1);
  end;

  pageind(pageind >= 4) = 4;
  %remember where the empty vectors are
  bad = pageind == 0;
  %but set the index so it won't give an error
  pageind(bad) = 1;

  %add 2 to skip over the pageind data itself
  pageind = (pageind-1)*nvec + 2;

  [i,j] = ndgrid(1:nx,1:ny);
  ind = sub2ind(size(data),i,j,pageind);
  u = data(ind);

  ind = sub2ind(size(data),i,j,pageind+1);
  v = data(ind);

  u(bad) = NaN;
  v(bad) = NaN;

  if (nvec == 3),
      ind = sub2ind(size(data),i,j,pageind+2);
      w = data(ind);
      w(bad) = NaN;
  end;
 case 2,                                % simple vector image
  u = A.Data(:,ry);
  v = A.Data(:,ry+ny);
 otherwise,
  % extract the correct vectors
  u = repmat(NaN,[nx ny]);
  v = repmat(NaN,[nx ny]);

  Achoice = A.Data(:,ry);               % defines where the correct
                                        % vectors are

  % don't ask me why this works this way
  for i = 1:6,
      if (i <= 4),
          off = (2*(i-1) + 1);
      else
          off = 7;
      end;

      mask = (Achoice == i);
      dat = A.Data(:,ry + off*ny);
      u(mask) = dat(mask);
      dat = A.Data(:,ry + (off+1)*ny);
      v(mask) = dat(mask);
  end;
end;

%flip and transpose everything around so that we end up with normal
%Matlab order for coordinates
y = y(end:-1:1);

%make normal plaid matrices
[x,y] = meshgrid(x,y);

u = u'*A.ScaleI(1) + A.ScaleI(2);
v = v'*A.ScaleI(1) + A.ScaleI(2);
u = flipud(u);
v = flipud(v);
if (nvec == 3),
    w = w'*A.ScaleI(1) + A.ScaleI(2);
    w = flipud(w);
end;

units.x = A.UnitX;
units.y = A.UnitY;
units.vel = A.UnitI;

% TODO: extract framing rate from M.Attributes

if (nvec == 3),
    varargout = {x,y,u,v,w};
else
    varargout = {x,y,u,v};
    if (nargout == 6),
        varargout{5} = [];
    end;
end;

if (nargout == length(varargout)+1),
    varargout{end+1} = units;
end;
