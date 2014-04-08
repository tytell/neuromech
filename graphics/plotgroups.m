function varargout = plotgroups(x,y,gp,mode,varargin)
%PLOTGROUPS    Plots data with multiple groupings
% plotgroups(x,y,gp,mode, etc...)
%   or plotgroups(y,gp,mode)
%   or plotgroups(x,y)
%
% Plots a data set that's catagorized by group.  Plots group membership in
% different ways.  Also can take means of data points all at the same x
% value.  x and y are the same size matrices containing the data points.
% x can also be a length m vector, and y is a m by n matrix.  The gp
% parameter specifies group membership and can contain multiple groups.
% It can either be a cell array of different group specifiers with the same
% size as y, or a 3D matrix where the first two dimensions are the same
% size as y, or, if y is a length n vector, it can be an n by p matrix.
% Mode specifies how different groups are identified.  It is a cell array
% of different specifications, which can contain
%
%   'm' - Vary the marker
%   'c' - Vary the color
%   'f' - Vary the fill
%
% Modes cannot repeat.
%
% For example, if
%    x = [1 2 3 4 1 2 3 1 5 6]';
%    y = [3 4 1 2 3 4 1 3 8 2]';
%    indiv = 'ABABABABAB';
%    speed = [1 1 1 1 2 2 2 2];
%
%    plotgroups(x,y,{indiv,speed},{'f','cm'});
%
% plots y vs x with individual coded by marker fill and speed coded by
% color and marker shape.
%
% Options:
%   'error' - Show error bars (can have values 'std','sem','all','none')
%   'means' - Take mean of points in the same group at the same x
%   'xOffset' - Define the offset of points with the same x, but different
%      group membership
%   'Markers' - Defines the list of markers for mode 'm'
%   'MarkerSize' - Defines the size of markers in points.
%   'Fill' - Defines the fill, or the start fill for mode 'f'
%   'Trace' - Connects points with the same group
%   'stdLineWidth' and 'semLineWidth' - Defines the width of lines for
%      'std' or 'sem' errors
%   'sortGroups' - Normally tries to keep the groups in the order they are
%      listed in the gp variable.  But if this is true, sort them
%      alphabetically or numerically
%   'legend' - Show a legend.  The normal legend function can't generally
%      handle plotgroups output very well.

% Mercurial revision hash: $Revision$ $Date$
% See https://bitbucket.org/tytell/matlab_variables/overview
% Copyright (c) 2011, Eric Tytell <tytell at jhu dot edu>

options.Error = 'sem';
options.xOffset = 'auto';
options.Markers = 'osd^v><ph';
options.MarkerSize = 8;
options.Fill = true;
options.stdLineWidth = 0.5;
options.semLineWidth = 0.5;
options.Color = [];
options.showTicks = [];
options.Trace = 'off';
options.isMarker = true;
options.means = false;
options.doMeans = logical([]);
options.sortGroups = false;
options.legend = false;
options.errorstyle = 'line';        % or 'bar'

optsynonyms = { ...
    {'xOffset',{'offx'}}, ...
    {'MarkerSize',{'ms'}}, ...
    {'stdLineWidth',{'stdlw'}}, ...
    {'semLineWidth',{'semlw'}}, ...
    {'showTicks',{'st'}}, ...
    {'sortGroups',{'sg'}}, ...
    };

if (nargin == 2)
    param = {x,y};
elseif (nargin == 3)
    param = {x,y,gp};
else
    param = {x,y,gp,mode};
end;

[options,param] = parsevarargin(options,[param varargin], 'leaveunknown', ...
    'synonyms',optsynonyms,'typecheck',false);

isX = true;
if (length(param) == 2),
    if (isnumeric(param{1}) && isnumeric(param{2})),
        x = param{1};
        y = param{2};
        gp = ones(size(y));
    elseif (iscell(param{2}) || ischar(param{2})),
        x = [];
        y = param{1};
        gp = param{2};
        isX = false;
    end;
    mode = {''};
elseif (length(param) == 3),
    if (iscell(param{2})),
        x = [];
        y = param{1};
        gp = param{2};
        mode = param{3};
        isX = false;
    else
        x = param{1};
        y = param{2};
        gp = param{3};
    end;
elseif (length(param) == 4),
    x = param{1};
    y = param{2};
    gp = param{3};
    mode = param{4};
end;

if (isX && ((ndims(x) ~= ndims(y)) || any(size(x) ~= size(y)))),
    error('x and y must be the same size');
end;
if (~iscell(gp)),
    if ((size(gp,1) ~= size(x,1)) || (size(gp,2) ~= size(x,2))),
        error('group must have the same first two sizes as x');
    end;
    if (length(mode) ~= size(gp,3)),
        error('mode must be as long as the number of pages in group.');
    end;
end;
if (iscell(gp)),
    for i = 1:length(gp),
        if ((ndims(gp{i}) ~= ndims(y)) || any(size(gp{i}) ~= size(y))),
            error('All elements of group must have the same size as y.');
        end;
    end;
    if (length(mode) ~= length(gp)),
        error('mode must be have as many entries as group.');
    end;
end;

% Process some of the options
if (options.doMeans),
    options.means = true;
end;
if (ischar(options.xOffset)),
    switch lower(options.xOffset),
     case 'auto',
      options.xOffset = [];
     case 'none',
      options.xOffset = 0;
    end;
end;
if (isempty(options.isMarker))
    if (ismember(options.Trace,{'on','meanonly'})),
        options.isMarker = false;
    else
        options.isMarker = true;
    end;
end;
if (ischar(options.Fill)),
    switch lower(options.Fill),
     case 'open',
      options.Fill = false;
     case 'filled',
      options.Fill = true;
    end;
end;

if (~isempty(x)),
    x0 = x;
    x = x(:);
end;

[m0 n0] = size(y);
sz = size(y);
y = y(:);
N = numel(y);

gpsz = size(gp);
if ((ndims(gp) == length(sz)+1) && all(gpsz(1:end-1) == sz(1:end-1))),
    gp = flatten(gp,1:ndims(gp)-1);
elseif ((ndims(gp) == length(sz)) && all(gpsz == sz)),
    gp = gp(:);
end;

if (isX),                       % check for repeated x values
    [xval,~,xind] = unique(x);
    nxval = length(xval);
    if (nxval == length(x)),
        isXRepeat = false;
    else
        isXRepeat = true;
    end;
else
    % deal with empty x
    x = zeros(size(y));
    xind = ones(size(y));
    nxval = 1;
    isXRepeat = true;
end;

good = find(isfinite(x) & isfinite(y));
xorig = x;
x = x(good);
y = y(good);
xind = xind(good);

% check the mode variable
mm = cat(2,mode{:});
if (length(unique(mm)) ~= length(mm)),
    error('Modes cannot repeat.');
end;

if (iscell(gp)),
    ngp = length(gp);
else
    ngp = size(gp,2);
end;

% by default, we don't do hsv colors, but we can, so set everything to
% default values
isHsv = 0;
hue = 0;
sat = 1;
val = 1;

% process the groups
for i = 1:ngp,
    if (iscell(gp)),
        gp1 = gp{i};
    else
        gp1 = gp(:,i);
    end;
    if (numel(gp1) ~= N),
        error(['Group specification variables must be the same size as ' ...
               'y.']);
    end;

    %find unique values
    %gpi1 will convert gp1 into group indices (1,2,...) instead of 
    %('A',B') or whatever
    %we don't want to remove any groups if all their y values are NaN.
    %So find unique group values in the entire matrix, then we'll remove
    %the NaN elements later
    [gpnames1,q,gpi1] = unique(gp1(:));
    gpi1 = shiftdim(gpi1);

    gpnum(i) = length(gpnames1);
    
    if (~options.sortGroups),
        %gpnames1 is sorted in alphabetical/numerical order.  We want it in
        %the order the values originally came in.
        %unfortunately sorting q does not give us the right order, so
        %we have to do something more complicated
        %look for the first time that each group index shows up in gpi1, then
        %sort those
        nameind = first(repmat(gpi1,[1 gpnum(i)]) == ...
            repmat(1:gpnum(i),[length(gpi1) 1]));
        [q,reshuffle] = sort(nameind);
        gpnames1 = gpnames1(reshuffle);
        %but now, the numbers in gpi1 don't correspond correctly to the
        %elements of gpnames1, so we need to reassign them
        %basically, we're reversing the sort on reshuffle and applying it to
        %gpi1
        [q,reshuffle] = sort(reshuffle);
        gpi1 = reshuffle(gpi1)';
    end;
    
    gpnames{i} = gpnames1;
    gpi(:,i) = gpi1;

    % build up the markers, colors, and fills from the mode specs
    for j = 1:length(mode{i}),
        switch mode{i}(j),
         case 'm',
          if (gpnum(i) > length(options.Markers)),
              options.Markers = repmat(options.Markers,...
                               [1 ceil(gpnum(i)/length(options.Markers))]);
          end;
         case 'f',
          isfill = [options.Fill 1-options.Fill];
          options.Fill = repmat(isfill,[1 ceil(gpnum(i)/2)]);
         case 'c',
          if (isempty(options.Color)),
              cmap = colormap;
              cind = round(linspace(1,size(cmap,1),gpnum(i)));
              options.Color = cmap(cind,:);
          else
              options.Color = shiftdim(options.Color);
              if (size(options.Color,1) < gpnum(i)),
                  options.Color = repmat(options.Color,...
                          [ceil(gpnum(i)/size(options.Color,1)) 1]);
              end;
          end;
         case 'h',
          hue = linspace(0,1,gpnum(i)+1);
          hue = hue(1:end-1);
          isHsv = 1;
         case 's',
          sat = linspace(0.5,1,gpnum(i));
          isHsv = 1;
         case 'v',
          val = linspace(0.5,1,gpnum(i));
          isHsv = 1;
        end;
    end;
end;

if (isempty(options.Color)),
    options.Color = [0 0 1];
end;

% Now find all the unique combinations of the different groups
% gpall identifies each individual combination
ngp = size(gpi,2);
[gpval,q,gpall] = unique(gpi,'rows');
gpi = gpi(good,:);
gpall = gpall(good);

ncomb = size(gpval,1);

%sort the group indices, so we can grab all the members of each group
%combination easily
[gpall,gpind] = sort(gpall);

gpi = gpi(gpind,:);
x = x(gpind);
xind = xind(gpind);
y = y(gpind);

% Set up the x offset
if (isXRepeat),
    if (isX),
        r = range(x(:));
    else
        r = 1;
    end;

    if (isempty(options.xOffset)),
        dx0 = r/(prod(gpnum)*nxval);
    
        incn = cumprod(gpnum);
        offxdist = [dx0 dx0*incn(1:end-1)];
    elseif (prod(size(options.xOffset)) == 1),
        dx0 = options.xOffset;
    
        incn = cumprod(gpnum);
        offxdist = [dx0 dx0*incn(1:end-1)];
    elseif (length(options.xOffset) == ngp),
        offxdist = options.xOffset;
    else
        error(['xOffset must be a scalar or a vector with one element per ' ...
               'group.']);
    end;

    if (ngp == 0),
        offx = zeros(size(x));
    else
        for i = 1:ngp,
            d = (1:gpnum(i))';
            d = d - (d(end)+d(1))/2;
            d = d*offxdist(i);
            
            offx(:,i) = d(gpi(:,i));
            offxval(:,i) = d(gpval(:,i));
        end;
    end;
    options.xOffset = offx;
end;

switch get(gca,'NextPlot'),
 case 'replace',
  cla reset;
  unhold = 1;
 case 'replacechildren',
  cla;
  unhold = 1;
 otherwise,
  unhold = 0;
end;

hold on;

% run through all the combinations
a = 1;
b = 1;
emptygroups = false(1,ncomb);

hStd = -1*ones(ncomb,1);
hSem = -1*ones(ncomb,1);

for i = 1:ncomb,
    while ((b <= length(gpall)) && (gpall(b) == i)),
        b = b+1;
    end;

    if (a <= length(gpall)),
        k = a:b-1;
    else
        k = [];
    end;
    x1 = x(k);
    if (~isempty(x1) & (size(options.xOffset,1) == length(x))),
        x1 = x1 + sum(options.xOffset(k,:),2);
    end;
    y1 = y(k);
    gpval1 = gpval(i,:);

    if (isempty(x1)),
        emptygroups(i) = true;
    end;

    isMean = options.means;
    if (isMean),
        [xind1,xord] = sort(xind(k));

        a2 = 1;
        b2 = 1;
        n = length(x1);

        xp = NaN([nxval 1]);
        yp = NaN([nxval 1]);
        stdp = NaN([nxval 1]);
        semp = NaN([nxval 1]);

        for j = 1:nxval,
            while ((b2 <= n) & (xind1(b2) == j)),
                b2 = b2+1;
            end;
            k2 = xord(a2:b2-1);

            if (isempty(k2)),
                xp(j) = NaN;
                yp(j) = NaN;
                stdp(j) = NaN;
                semp(j) = NaN;
            else
                xp(j) = x1(xord(a2));

                y2 = y1(k2);
                yp(j) = mean(y2);
                stdp(j) = std(y2);
                semp(j) = stdp(j)/sqrt(length(y2));
            end;
            a2 = b2;
        end;

        q = find(isfinite(yp) & isfinite(xp));
        x1 = xp(q);
        y1 = yp(q);
        std1 = stdp(q);
        sem1 = semp(q);
    end;

    col1 = options.Color(1,:);
    marker1 = options.Markers(1);
    fill1 = options.Fill(1);
    if (isHsv),
        hue1 = hue(1);
        sat1 = sat(1);
        val1 = val(1);
    end;
    for j = 1:ngp,
        for k = 1:length(mode{j}),
            switch mode{j}(k),
             case 'c',
              col1 = options.Color(gpval1(j),:);
             case 'm',
              marker1 = options.Markers(gpval1(j));
             case 'f',
              fill1 = options.Fill(gpval1(j));
             case 'h',
              hue1 = hue(gpval1(j));
             case 's',
              sat1 = sat(gpval1(j));
             case 'v',
              val1 = val(gpval1(j));
            end;
        end;
    end;

    if (isHsv),
        col1 = hsv2rgb(hue1,sat1,val1);
    end;

    if (isMean),
        if (strcmpi(options.Error,'std') || strcmpi(options.Error,'all')),
            hStd(i) = -1;
            if (ismember(options.Trace,{'on'}) && strcmpi(options.errorstyle,'bar')),
                px = [x1; x1(end:-1:1)];
                py = [y1; y1(end:-1:1)] + [std1; -std1(end:-1:1)];
                h1 = fill(px,py,col1,'EdgeColor','none');
                hStd(i) = h1;
            else
                px = [x1 x1 NaN([length(x1) 1])]';
                py = [y1 y1 NaN([length(y1) 1])]' + ...
                     [std1 -std1 NaN([length(y1) 1])]';
            end;
            if (~isempty(px) && ~ishandle(hStd(i))),
                hStd(i) = plot(px(:),py(:),'Color',col1,...
                               'LineWidth',options.stdLineWidth);
            elseif ~ishandle(hStd(i))
                %have to use line here to force the generation of a
                %graphics handle even for a trace with no points
                hStd(i) = line('XData',[],'YData',[],'Color',col1,...
                               'LineStyle','none');
            end;
        end;
        if (strcmpi(options.Error,'sem') || strcmpi(options.Error,'all')),
            if (ismember(options.Trace,{'on'})),
                px = [x1; x1(end:-1:1)];
                py = [y1; y1(end:-1:1)] + [sem1; -sem1(end:-1:1)];
                h1 = fill(px,py,col1,'EdgeColor','none');
                hSem(i) = h1;
            else
                px = [x1 x1 NaN([length(x1) 1])]';
                py = [y1 y1 NaN([length(y1) 1])]' + ...
                     [sem1 -sem1 NaN([length(y1) 1])]';
            end;
            if (~isempty(px) && ~ishandle(hSem(i))),
                hSem(i) = plot(px(:),py(:),'Color',col1,...
                               'LineWidth',options.semLineWidth);
            elseif ~ishandle(hSem(i))
                hSem(i) = line('XData',[],'YData',[],'Color',col1,...
                               'LineStyle','none');
            end;
        end;

    end;

    if (isempty(options.isMarker) || ~options.isMarker),
        marker1 = '';
    end;
    if (ismember(options.Trace,{'on','meanonly'})),
        marker1 = ['-' marker1];
    end;

    if (~isempty(marker1)),
        if (fill1),
            markeropts = {'MarkerFaceColor',col1,...
                          'MarkerEdgeColor','none','Color',col1, ...
                          'MarkerSize',options.MarkerSize};
        else
            markeropts = {'MarkerFaceColor','none',...
                          'MarkerEdgeColor',col1,'Color',col1, ...
                          'MarkerSize',options.MarkerSize};
        end;

        if (~isempty(x1)),
            hpt(i) = plot(x1,y1,marker1,markeropts{:});
        else
            hpt(i) = line('XData',x1,'YData',y1,'Marker','none', ...
                          'LineStyle','none');
        end;
    end;

    a = b;
end;

if (~isX),
    %stn contains the indices of the groups we want to generate ticks for
    if (isempty(options.showTicks)),
        stn = 1:ngp;
    else
        stn = options.showTicks;
    end;

    %get the unique combinations of the selected groups
    gptick = unique(gpval(:,stn),'rows');

    %build up the names of the combinations
    for i = 1:length(gptick),
        xtickName = '';

        %run through all the groups we're showing ticks for
        for mm = 1:length(stn),
            m = stn(mm);

            %get the name of this group (number in gptick)
            switch class(gpnames{m}),
             case 'char',
              xtn = gpnames{m}(gptick(i,m));
             case 'cell',
              xtn = gpnames{m}{gptick(i,m)};
             otherwise,
              xtn = num2str(gpnames{m}(gptick(i,m)));
            end;
            %add it to the end of xtickName
            xtickName = [xtickName xtn ' '];
        end;

        %remove the final space
        xtickName = xtickName(1:end-1);
        xtickNames{i} = xtickName;

        %find all x values that have this combination of groups
        k = find(all(gpval(:,stn) == repmat(gptick(i,:),[ncomb 1]), 2));
        %and get the mean position - we need to do this because we might
        %be showing the name of one group centered under several
        %subgroups
        xticks(i) = mean(sum(offxval(k,:),2));
    end;

    %make sure we haven't produced the same tick position multiple times
    [xticks,ind] = unique(xticks);
    xtickNames = xtickNames(ind);

    set(gca, 'XTick',xticks,'XTickLabel',xtickNames);
end;

combnames = cell(ncomb,1);
for i = 1:ncomb,
    if (ngp > 1),
        combname1 = '(';
    else
        combname1 = '';
    end;
    for j = 1:ngp,
        switch class(gpnames{j}),
          case 'char',
            combname1 = strcat(combname1,gpnames{j}(gpval(i,j)),',');
          case 'cell',
            combname1 = strcat(combname1,gpnames{j}{gpval(i,j)},',');
          otherwise,
            combname1 = strcat(combname1,num2str(gpnames{j}(gpval(i,j))),...
                             ',');
        end;
    end;
    if (ngp > 1),
        combname1(end) = ')';
    else
        combname1 = combname1(1:end-1);
    end;
    combnames{i} = combname1;
end;

if (any(emptygroups)),
    fprintf('%d empty groups found:\n  ',sum(emptygroups));
    fprintf('%s ',combnames{emptygroups});
    fprintf('\n');
end;

if (options.legend),
    legend(hpt(~emptygroups),combnames(~emptygroups));
end;

if (unhold),
    hold off;
end;
        
if (nargout == 1),
    varargout{1} = hpt;
end;
