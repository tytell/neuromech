function varargout = xyzvec(x,y,z,varargin)

switch nargin,
    case 1,
        varargin = {x};
    case 2,
        varargin = {x,y};
    case 3,
        varargin = {x,y,z};
    otherwise
        varargin = [x y z varargin];
end;

sz(1,:) = cellfun(@(x) size(x,1), varargin);
sz(2,:) = cellfun(@(x) size(x,2), varargin);
sz(3,:) = cellfun(@(x) size(x,3), varargin);

if ((nargin >= 2) && isnumeric(varargin{1}) && ...
        ((sz(1,1) == 2) || (sz(1,1) == 3)) && ...
        isnumeric(varargin{2}) && (ndims(varargin{2}) == 2) && ...
        (sum(sz(:,2) > 1) == 1))
    [X,sz] = deal(varargin{1:2});
    tovec = false;
    p = 3;
elseif ((nargin >= 3) && all(cellfun(@isnumeric,varargin(1:3))) && ...
        all(all(sz(:,1:3) == sz(:,[1 1 1]))))
    [x,y,z] = deal(varargin{1:3});
    tovec = true;
    p = 4;
elseif ((nargin >= 2) && all(cellfun(@isnumeric,varargin(1:2))) && ...
        all(all(sz(:,1:2) == sz(:,[1 1]))))
    [x,y] = deal(varargin{1:2});
    z = [];
    tovec = true;
    p = 3;
elseif ((nargin >= 1) && isnumeric(varargin{1}))
    X = varargin{1};
    sz = size(X);
    if (length(sz) == 2)
        sz = [sz(2) 1];
    else
        sz = sz(2:end);
    end;
    tovec = false;
    p = 2;
end;

% opt = parsevarargin(opt,varargin(p:end),p);

if (tovec)
    sz = size(x);
    n = numel(x);
    
    x = reshape(x,[1 n]);
    y = reshape(y,[1 n]);
    if (isempty(z))
        X = cat(1,x,y);
    else
        z = reshape(z,[1 n]);
        X = cat(1,x,y,z);
    end;
    varargout = {X, sz};
else
    Xp = cell(1,size(X,1));
    for i = 1:size(X,1),
        Xp{i} = squeeze(X(i,:,:));
        Xp{i} = reshape(Xp{i},sz);
    end;
    varargout = Xp;
end;
