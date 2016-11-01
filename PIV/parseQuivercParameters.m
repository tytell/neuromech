function [x,y,z,u,v,w,col,opts] = parseQuivercParameters(param)
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell

% get the size of all the arguments
sz(:,1) = cellfun('size',param,1)';
sz(:,2) = cellfun('size',param,2)';
sz(:,3) = cellfun('size',param,3)';

% assume we're 2D at first
opts.is3D = 0;
z = [];
w = [];

% check for various combinations of x,y,u,v
if ((length(param) >= 6) && all(cellfun('isclass',param(1:6),'double')) && ...
    all(sz(1,:) == sz(2,:)) && ...
    all(sz(1,:) == sz(3,:)) && all(sz(1,:) == sz(4,:)) && ...
    all(sz(1,:) == sz(4,:)) && all(sz(1,:) == sz(5,:))),
    x = param{1};
    y = param{2};
    z = param{3};
    u = param{4};
    v = param{5};
    w = param{6};
    opts.is3D = 1;
    p = 7;
elseif ((length(param) >= 6) && ...
        all(cellfun('isclass',param(1:6),'double')) && ...
        all(sz(4,:) == sz(5,:)) && ...
        (prod(sz(1,:)) == sz(4,2)) && (prod(sz(2,:)) == sz(4,1)) && ...
        (prod(sz(3,:)) == sz(4,3))),
    % quiverc(x(1,:,1),y(:,1,1),z(1,1,:),u,v,w,...)
    [x,y,z] = meshgrid(param{1},param{2},param{3});
    u = param{4};
    v = param{5};
    w = param{6};
    opts.is3D = 1;
    p = 7;
elseif ((length(param) >= 4) && ...
        all(cellfun('isclass',param(1:4),'double')) && ...
        all(sz(1,:) == sz(2,:)) && ...
        all(sz(1,:) == sz(3,:)) && all(sz(1,:) == sz(4,:))),
% quiverc(x,y,u,v,...)
    x = param{1};
    y = param{2};
    u = param{3};
    v = param{4};
    p = 5;
elseif ((length(param) >= 4) && ...
        all(cellfun('isclass',param(1:4),'double')) && ...
        all(sz(3,:) == sz(4,:)) && ...
        (prod(sz(1,:)) == sz(3,2)) && (prod(sz(2,:)) == sz(3,1))),
% quiverc(x(1,:),y(:,1),u,v,...)
    [x,y] = meshgrid(param{1},param{2});
    u = param{3};
    v = param{4};
    p = 5;
elseif ((length(param) >= 2) && ...
        all(cellfun('isclass',param(1:2),'double')) && ...
        all(sz(1,:) == sz(2,:))),
% quiverc(u,v,...)
    u = param{1};
    v = param{2};
    [x,y] = meshgrid(1:size(u,2),1:size(u,1));
    
    p = 3;
else
    error('Parameter sizes do not match');
end;

if (~opts.is3D),
    if (any(sz(:,3) > 1)),
        error('Cannot handle 3D matrices in 2D plots.');
    end;
    sz = sz(:,1:2);                     % remove third dimension
end;

opts.Color = [];
if ((p <= length(param)) && ~ischar(param{p}) && ...
    (ndims(param{p}) == ndims(u)) && (all(size(param{p}) == size(u)))),
    opts.Color = param{p};
    p = p+1;
    opts.isColMatrix = 1;
elseif ((p <= length(param)) && ~ischar(param{p}) && ...
        (numel(param{p}) == 3)),
    opts.Color = param{p};
    p = p+1;
    opts.isColMatrix = 0;
else
    opts.isColMatrix = 0;
end;

% process all the options
param = param(p:end);

[opts,param] = matchQuivercOption(opts,param,'ScaleAsPrevious','ps');
if (opts.ScaleAsPrevious),
    obj = findobj(gca,'Tag','quiverc');
    if (~isempty(obj)),
        opts0 = get(obj(1),'UserData');
        opts.Scale = opts0.Scale;
        opts.HeadRange = opts0.HeadRange;
        opts.AbsHeadRange = opts0.AbsHeadRange;
        opts.ScaleAsPrevious = 1;
    end;
end;

[opts,param,m1] = matchQuivercOption(opts,param,'ScaleFactor','s',1);
[opts,param,m2] = matchQuivercOption(opts,param,'RelScale','rs',0.95);
[opts,param,m3] = matchQuivercOption(opts,param,'AbsScale','as',[]);
if (m3 && (m1 || m2)),
    warning('quiverc:AbsScaleOverride','AbsScale option overrides ScaleFactor or RelScale option.');
end;

[opts,param,m1] = matchQuivercOption(opts,param,'ScaleRange','sr',[0 1]);
[opts,param,m2] = matchQuivercOption(opts,param,'AbsScaleRange',...
                                     'asr',[]);
if (m1 || m2),
    isLenScaleRange = 1;
    if (m1 && m2),
        warning('quiverc:AbsScaleOverride','AbsScaleRange option overrides ScaleRange option.');
    end;
else
    isLenScaleRange = 0;
end;

[opts,param] = matchQuivercOption(opts,param,'HeadSize','hs',0.4);
[opts,param,m1] = matchQuivercOption(opts,param,'HeadRange','hr',[0 1]);
[opts,param,m2] = matchQuivercOption(opts,param,...
                                     'AbsHeadRange','ahr',[]);
if (~m1 && ~m2 && isLenScaleRange),
    opts.HeadRange = opts.ScaleRange;
elseif (m1 && m2),
    warning('quiverc:AbsHeadRangeOverride','AbsHeadRange option overrides HeadRange option.');
end;

[opts,param,m1] = matchQuivercOption(opts,param,'Truncate','t',1);
if (m1),
    warning('quiverc:truncateIsObsolete',...
        'Obsolete option ''truncate''.  Use ''ScaleRange'' instead.');
    opts = rmfield(opts,'Truncate');
    opts.ScaleRange = [0 trunc];
end;
[opts,param] = matchQuivercOption(opts,param,'Show','sh',1);

[opts,param] = matchQuivercOption(opts,param,'NoHeads','nh');
[opts,param] = matchQuivercOption(opts,param,'NoTails','nt');
[opts,param] = matchQuivercOption(opts,param,'Angle','ang');
if (opts.Angle),
    opts.NoHeads = 1;
end;
if (opts.NoHeads && opts.NoTails),
    error('You must either display heads or tails.');
end;
[opts,param] = matchQuivercOption(opts,param,'HeadFacets3d','hf3d',8);
[opts,param] = matchQuivercOption(opts,param,'LineWidth','lw',0.5);

[opts,param] = matchQuivercOption(opts,param,'CorrectAspectRatio','car');
[opts,param] = matchQuivercOption(opts,param,'Polar','pol');

if exist('gobjects', 'file') == 2
    axdef = gobjects(0);
else
    axdef = [];
end
[opts,param] = matchQuivercOption(opts,param,'Axes','ax',axdef);

if (~opts.isColMatrix),
    [opts,param,m] = matchQuivercOption(opts,param,'Color','','');
    if (strcmpi(opts.Color,'cycle')),
        opts.Color = 'b';
        opts.isColMatrix = 0;

        obj = findobj(gca,'Tag','quiverc');
        if (~isempty(obj)),
            switch get(obj(1),'Type'),
             case 'line',
              prevcol = get(obj(1),'Color');
             case 'patch',
              ecol = get(obj(1),'EdgeColor');
              if (isnumeric(ecol)),
                  prevcol = ecol;
              else
                  prevcol = get(obj(1),'FaceColor');
              end;
            end;
            colororder = get(gca,'ColorOrder');

            if (all(size(prevcol) == [1 3])),
                ind = find(all(colororder == ...
                               repmat(prevcol,[size(colororder,1) 1]),2));
                if (~isempty(ind)),
                    ind = mod(ind+1,size(colororder,1));
                    col = colororder(ind,:);
                    opts.Color = col;
                end;
            end;
        end;
    end;
end;
    
colornames = {'y','yellow'; 'm','magenta'; 'c','cyan'; 'r','red'; ...
              'g','green'; 'b','blue'; 'w','white'; 'k','black'};

% Check for color matrix or color string
if (~isempty(param) && ischar(param{1}) && ...
    (numel(strmatch(param{1},colornames(:),'exact')) == 1)),
    if (~isempty(opts.Color)),
        warning('quiverc:ColorOverride','Color option overrides color variable.');
    end;
    opts.Color = param{1};
    param = param(2:end);
    opts.isColMatrix = 0;
    p = p+1;
end;
if (isempty(opts.Color)),
    opts.Color = 'b';
    opts.isColMatrix = 0;
end;
col = opts.Color;

% error if there's anything left
if (~isempty(param)),
  isc = cellfun('isclass',param,'char');

  str = sprintf('%s, ',param{isc});
  str = str(1:end-2);

  if (sum(isc) == length(param)),
      error('Unrecognized options %s', str);
  elseif (isempty(isc)),
      error('%d unrecognized parameters', length(param));
  else
      error(['%d unrecognized parameters and %d unrecognized ' ...
                     'options (%s)'], length(param)-length(isc), ...
                    length(isc), str);
  end
end;

opts.x = x;
opts.y = y;
opts.z = z;
opts.u = u;
opts.v = v;
opts.w = w;
