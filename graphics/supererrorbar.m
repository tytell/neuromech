function [h,err] = supererrorbar(varargin)
% function [h,err] = supererrorbar(x,y,gp)
%         or         supererrorbar(x,y)
%         or         supererrorbar(y)
%
% plus options
%    supererrorbar(...,linespec)
%       where linespec is the standard linespec from plot
%    'sem','std',('rng' or 'range'): plot the corresponding thing as
%        the error

k = find(cellfun('isclass',varargin,'char'));

if (~isempty(k)),
    opt = k(1);
else,
    opt = nargin+1;
end;

%check number of parameters before character options
switch opt-1,
 case 1,
  y = varargin{1};
  x = repmat(1:size(y,2),[size(y,1) 1]);
  gp = ones(size(y));
 case 2,
  [x,y] = deal(varargin{1:2});
  gp = ones(size(y));
 case 3,
  [x,y,gp] = deal(varargin{1:3});
 otherwise,
  error('Unknown parameters');
end;

if (any(size(x) == 1)),
    if (length(x) ~= size(y,2)),
        error('Vector x must have the same number of columns as y.');
    end;
    x = repmat(reshape(x,[1 size(y,2)]),[size(y,1) 1]);
end;

if ((ndims(x) ~= ndims(y)) | (ndims(gp) ~= ndims(y)) | ...
    any(size(x) ~= size(y)) | any(size(gp) ~= size(y))),
    error('x,y,and gp matrices must be the same size');
end;

%set up the defaults
options = struct('Color',{''},'Marker',{'o'},'LineStyle',{'-'},...
              'Error',{@sebSEM});

%check if the first option has the same size as y.  If so, it's a group,
%not an option
if ((opt <= nargin) & (ndims(varargin{opt}) == ndims(y)) & ...
    all(size(varargin{opt}) == size(y))),
    gp = varargin{opt};
    opt = opt+1;
end;

%process a full LineSpec, if necessary
if ((opt <= nargin) & ischar(varargin{opt}) & ...
    (length(varargin{opt}) <= 4)),
    col = regexp(varargin{opt},'[ymcrbwk]');
    mark = regexp(varargin{opt},'[+o*.xsd^v<>ph]');
    [ln1,ln2] = regexp(varargin{opt},'(--|-\.|-|:)','once');

    %check whether there's anything left over, after we've matched all of
    %these
    match = unique([col mark ln1 ln2]);
    if (length(match) == length(varargin{opt})),
        %it's OK
        if (~isempty(col)),
            options.Color = varargin{opt}(col);
        end;
        if (~isempty(mark)),
            options.Marker = varargin{opt}(mark);
        end;
        if (~isempty(ln1)),
            options.LineStyle = varargin{opt}(ln1:ln2);
        end;
        opt = opt+1;
    end;
end;

if (opt <= nargin),
    opts = varargin(opt:end);

    %and process the remaining options
    i = 1;
    while (i <= length(opts)),
        switch lower(opts{i}),
         case 'color',
          options.Color = opts{i+1};
          i = i+2;
         case 'sem',
          options.Error = @sebSEM;
          i = i+1;
         case 'std',
          options.Error = @std;
          i = i+1;
         case {'rng','range'},
          options.Error = @range;
          i = i+1;
         otherwise,
          error('Unrecognized option %s.',opts{i});
        end;
    end;
end;

if (isnumeric(gp)),
    nonan = find(isfinite(x) & isfinite(y) & isfinite(gp));
else
    nonan = find(isfinite(x) & isfinite(y));
end;
x = x(nonan);
y = y(nonan);
gp = gp(nonan);

%get unique groups
[gpval,q,gpind] = unique(gp(:));
ngp = length(gpval);

%and unique positions
[xval,q,xind] = unique(x(:));
nx = length(xval);
xval = repmat(xval,[1 ngp]);

%offsets for different groups in a single x
xoff = range(xval)/nx/ngp;

for j = 1:ngp,
    for i = 1:nx,
        k = find((xind == i) & (gpind == j));
        mn(i,j) = mean(y(k));
        err(i,j) = feval(options.Error,y(k));
    end;
end;

errx = permute(repmat(xval,[1 1 2]),[3 1 2]);
errx(3,:,:) = NaN;
errx = reshape(errx,[3*size(xval,1) size(xval,2)]);
erry = permute(repmat(mn,[1 1 2]),[3 1 2]);
erry(1,:,:) = erry(1,:,:) - permute(err,[3 1 2]);
erry(2,:,:) = erry(2,:,:) + permute(err,[3 1 2]);
erry(3,:,:) = NaN;
erry = reshape(erry,[3*size(xval,1) size(xval,2)]);

%check whether the plot is held
isheld = ishold;

if (~isempty(options.Color)),
    colopt = {'Color',options.Color};
else
    colopt = {};
end;
h = plot(xval,mn,colopt{:},'Marker',options.Marker,'LineStyle','none');
hold on;
h = cat(1,h,plot(errx,erry,colopt{:},'LineStyle',options.LineStyle));

if (~isheld),
    hold off;
end;



%--------------------------------------------------
function sem = sebSEM(y)

sem = std(y(:))./sqrt(numel(y));




