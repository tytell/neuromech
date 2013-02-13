function w = vorticity(varargin)
% function w = vorticity(x,y,u,v,err)

if (nargin == 1),
    %passed in a DU parameter from velderiv
    dudy = varargin{1}.dudy;
    dvdx = varargin{1}.dvdx;
else,
    x = varargin{1};
    y = varargin{2};
    u = varargin{3};
    v = varargin{4};

    if (nargin == 5)
        err = varargin{5};

        k = find(err == -1);
        u(k) = NaN;
        v(k) = NaN;
    end;

    if ((size(x,3) == 1) & (size(u,3) > 1)),
        x = repmat(x,[1 1 size(u,3)]);
        y = repmat(y,[1 1 size(u,3)]);
    end;
    
    dudy = repmat(NaN,size(u));
    dvdx = repmat(NaN,size(v));

    dudy(2:end-1,:,:) = (u(3:end,:,:) - u(1:end-2,:,:)) ./ ...
        (y(3:end,:,:) - y(1:end-2,:,:));
    dvdx(:,2:end-1,:) = (v(:,3:end,:) - v(:,1:end-2,:)) ./ ...
        (x(:,3:end,:) - x(:,1:end-2,:));
end;

% see Faber 1995 p. 22 for sign
w = dvdx - dudy;