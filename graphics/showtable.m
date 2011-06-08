function h = showtable(varargin)
% h = showtable(tab, ...)
%   or  showtable(ax,tab, ...)
%
% The main parameter is a matrix (numeric or cell) that will be displayed
% in a table.  The table can contain numbers, strings, or small plots (more
% or less like "sparklines" as described by Tufte
% see http://www.edwardtufte.com/bboard/q-and-a-fetch-msg?msg_id=0001OR).
%
% There are many options for adjusting fonts, alignments, number formats,
% and spacing.  Most options can be specified as a single value (for all
% cells), as a row specifying different values for each column, as a column
% specifying different values for each row, or as a full matrix with
% different values for each cell.
%
% Options:
%    'mode' -
%       'plot' - Plots the table in Matlab's figure
%       'latex' - Attempts to produce LaTeX code for the table
%    'rownames', 'colnames' - Cell strings with names for the rows or
%      columns
%    'font' - Font properties for the main text.  Can be just a font name
%      or a cell array of property names and values that the text function
%      can handle.
%    'colnamesfont','rownamesfont' - Font properties for the column names
%      or row names.  Can be just a font name or a cell array with property
%      names and values.
%    'align' - Describes the alignment within cells with two characters.
%      For horizontal alignment, values can be 'l', 'c', or 'r', for left,
%      center, or right.  For vertical alignment, values can be 't','m', or
%      'b' for top, middle, or bottom.  Default is 'bc' for aligned to bottom
%      and centered horizontally.  Values can come in any order.
%    'colwidthmode','rowheightmode' - Can be
%      'tight' - minimum width or height
%      'even' - same width for all in which the largest item fits
%      'manual' - widths or heights specified in 'colwidth' or
%                'rowheight' in points
%      'span' - cover the entire axes.  'colwidth' or 'rowheight'
%               parameters then specify the relative size of each column
%               (e.g., showtable(tab,'colwidthmode','span','colwidth',[1 1
%               1 2]) would give the last column twice the size of the
%               others).  Any columns given zero size expand to fill the
%               available space.  plot cells by default have zero size.
%    'colwidth','rowheight' - Specifies the width or height of the columns
%      or rows in points (for 'manual' mode) or relative sizes (for 'span'
%      mode)
%    'rowgap','colgap' - Size of the gap between rows or columns in points.
%    'xlim','ylim' - X or Y limits for plots.  Can be
%      'linkrows','linkcols' - Link the limits along columns or rows
%      'tight' - Tight limits individually for cells
%      2 element matrix - Manual specifications for all cells
%      cell array of 2 element matrices - Manual specifications for
%        columns, rows, or each cell individually
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell <tytell at jhu dot edu>

opt.rownames = {};
opt.colnames = {};
opt.font = 'Times New Roman';
opt.colnamesfont = {'FontAngle','italic'};
opt.colnamesheight = [];
opt.rownamesfont = {'FontWeight','bold'};
opt.align = 'bc';
opt.format = {};
opt.colwidthmode = 'even';
opt.rowheightmode = 'even';
opt.colwidth = [];
opt.rowheight = [];
opt.rowgap = 0;
opt.colgap = 0;
opt.xlim = 'tight';
opt.ylim = 'tight';
opt.mode = 'plot';

%first look for axes options
origax = gca;
if ((nargin >= 2) && (numel(varargin{1}) == 1) && ...
        ishandle(varargin{1}) && strcmp(get(varargin{1},'Type'),'axes')),
    [ax,X] = varargin{1:2};
    p = 3;
else
    ax = gca;
    X = varargin{1};
    p = 2;
end;

if (~strcmp(get(ax,'NextPlot'),'add'))
    cla(ax,'reset');
end;

%parse the options
opt = parsevarargin(opt, varargin(p:end), p, 'typecheck',false,...
    'multival',{'mode',{'plot','latex'}},'exact');

%set up the X matrix
if (ndims(X) ~= 2),
    error('showtable:argsize','X must be two dimensional');
end;
if (isnumeric(X)),
    X = num2cell(X);
end;

nr = size(X,1);
nc = size(X,2);

%start putting together the texprops matrix
textprops = cell(size(X));
if (ischar(opt.font)),
    %single value -> put it everywhere
    [textprops{:}] = deal({'FontName',opt.font});
elseif (iscellstr(opt.font)),
    [textprops{:}] = deal(opt.font);
elseif (iscell(opt.font)),
    if ((size(opt.font,1) == 1) && (length(opt.font) == size(X,2))),
        textprops = repmat(opt.font(:)',[size(X,1) 1]);
    elseif ((size(opt.font,2) == 1) && (length(opt.font) == size(X,1))),
        textprops = repmat(opt.font(:),[1 size(X,2)]);
    elseif ((ndims(opt.font) == 2) && (all(size(opt.font) == size(X)))),
        textprops = opt.font;
    end;
end;

%format matrix
numfmt = '%g';
if (isempty(opt.format)),
    fmt = cell(size(X));
elseif (ischar(opt.format)),
    fmt = cell(size(X));
    numfmt = opt.format;
elseif (iscell(opt.format) && (length(opt.format) == nr)),
    fmt = repmat(opt.format(:),[1 nc]);
elseif (iscell(opt.format) && (length(opt.format) == nc)),
    fmt = repmat(opt.format(:)',[nr 1]);
elseif (iscell(opt.format) && (ndims(opt.format) == 2) && ...
        all(size(opt.format) == size(X))),
    fmt = opt.format;
else
    error('showtable:argsize','format matrix is the wrong size');
end;

for i = 1:nr,
    for j = 1:nc,
        if (ischar(X{i,j})),
            if (isempty(fmt{i,j}))
                fmt{i,j} = '%s';
            end;
        elseif (isempty(X{i,j}))
            fmt{i,j} = '';
        elseif ((isnumeric(X{i,j}) || islogical(X{i,j})) && ...
                (numel(X{i,j}) == 1)),
            if (isempty(fmt{i,j}))
                fmt{i,j} = numfmt;
            end;
        elseif (isnumeric(X{i,j})),
            if (isempty(fmt{i,j}))
                fmt{i,j} = @plot;
            end;
            X{i,j} = {X{i,j}};
        elseif (iscell(X{i,j}) && isa(X{i,j}{1},'function_handle')),
            fmt{i,j} = X{i,j}{1};
            X{i,j} = X{i,j}(2:end);
        elseif (iscell(X{i,j})),
            if (isempty(fmt{i,j}))
                fmt{i,j} = @plot;
            end
        end;
    end;
end;

%x limits
if (ischar(opt.xlim)),
    xlim = cell(nr,nc);
    [xlim{:}] = deal(opt.xlim);
elseif (iscell(opt.xlim) && (length(opt.xlim) == nc)),
    xlim = repmat(opt.xlim(:)',[nr 1]);
elseif (iscell(opt.xlim) && (ndims(opt.xlim) == 2) && ...
        all(size(opt.xlim) == size(X))),
    xlim = opt.xlim;
else
    error('Unrecognized xlim option');
end;

%y limits
if (ischar(opt.ylim)),
    ylim = cell(nr,nc);
    [ylim{:}] = deal(opt.ylim);
elseif (iscell(opt.ylim) && (length(opt.ylim) == nr)),
    ylim = repmat(opt.ylim(:),[1 nc]);
elseif (iscell(opt.ylim) && (ndims(opt.ylim) == 2) && ...
        all(size(opt.ylim) == size(X))),
    ylim = opt.ylim;
else
    error('Unrecognized xlim option');
end;

%row names
if (~isempty(opt.rownames)),
    if (length(opt.rownames) ~= nr),
        error('showtable:argsize','rownames must have the same number of rows as X');
    end;
    X = cat(2,opt.rownames(:),X);
    
    [ok, opt.rownamesfont] = parseheaderfont(nr,opt.rownamesfont, textprops{1,1});
    if (~ok),
        error('showtable:rownamesfont','Cannot understand rownamesfont');
    end;
    textprops = cat(2,opt.rownamesfont',textprops);
    
    fmt = cat(2,cell(nr,1),fmt);
    [fmt{:,1}] = deal('%s');
    
    xlim = cat(2,cell(nr,1),xlim);
    ylim = cat(2,cell(nr,1),ylim);
    
    isrownames = true;
    
    a = 2;
else
    isrownames = false;
    a = 1;
end;

%column names
if (~isempty(opt.colnames))
    if (length(opt.colnames) ~= nc),
        error('showtable:argsize','colnames must have the same number of columns as X');
    end;
    X = cat(1,cell(1,size(X,2)),X);
    X(1,a:end) = opt.colnames;

    fmt = cat(1,cell(1,size(fmt,2)),fmt);
    [fmt{1,a:end}] = deal('%s');
    
    [ok, opt.colnamesfont] = parseheaderfont(nc,opt.colnamesfont, textprops{1,a});
    if (~ok),
        error('showtable:colnamesfont','Cannot understand colnamesfont');
    end;
    textprops = cat(1,cell(1,size(textprops,2)), textprops);
    textprops(1,a:end) = opt.colnamesfont;
    
    xlim = cat(1,cell(1,size(xlim,2)),xlim);
    ylim = cat(1,cell(1,size(ylim,2)),ylim);
    
    iscolnames = true;
else
    iscolnames = false;
end;

nr = size(X,1);
nc = size(X,2);

%alignment
if (ischar(opt.align)),
    opt.align = repmat({opt.align},[nr nc]);
elseif (iscellstr(opt.align) && (length(opt.align) == nc)),
    opt.align = repmat(opt.align(:)',[nr 1]);
elseif (iscellstr(opt.align) && (length(opt.align) == nr)),
    opt.align = repmat(opt.align(:),[1 nc]);
elseif (iscell(opt.align) && (ndims(opt.align) == 2) && ...
        all(size(opt.align) == size(X))),
    %everything's OK
else
    error('showtable:argsize','Wrong size alignment matrix');
end;

offx = zeros(size(X));
offy = zeros(size(X));
for i = 1:nr,
    for j = 1:nc,
        if (length(opt.align{i,j}) ~= 2),
            error('showtable:argsize','Alignment matrix error at position %d,%d',i,j);
        end;
        ishoriz = ismember(lower(opt.align{i,j}),'lcr');
        isvert = ismember(lower(opt.align{i,j}),'tmb');
        
        switch lower(opt.align{i,j}(ishoriz)),
            case 'l',
                textprops{i,j} = [textprops{i,j} {'HorizontalAlignment','left'}];
            case 'c',
                offx(i,j) = 0.5;
                textprops{i,j} = [textprops{i,j} {'HorizontalAlignment','center'}];
            case 'r',
                offx(i,j) = 1;
                textprops{i,j} = [textprops{i,j} {'HorizontalAlignment','right'}];
            otherwise,
                error('Unrecognized aligment option %s at position %d,%d',...
                    opt.align{i,j}(ishoriz), i,j);
        end;                
        switch lower(opt.align{i,j}(isvert)),
            case 't',
                offy(i,j) = 1;
                textprops{i,j} = [textprops{i,j} {'VerticalAlignment','top'}];
            case 'm',
                offy(i,j) = 0.5;
                textprops{i,j} = [textprops{i,j} {'VerticalAlignment','middle'}];
            case 'b',
                textprops{i,j} = [textprops{i,j} {'VerticalAlignment','bottom'}];
            otherwise,
                error('Unrecognized aligment option %s at position %d,%d',...
                    opt.align{i,j}(isvert), i,j);
        end;                
    end;
end;

%margin around the cells
if (numel(opt.rowgap) == 1),
    rowgap = opt.rowgap * ones(nr+1,1);
    rowgap([1 nr+1]) = 0;
elseif (length(opt.rowgap) == nr-1),
    rowgap = [0; opt.rowgap(:); 0];
elseif (length(opt.rowgap) == nr),
    rowgap = [opt.rowgap(:); 0];
else
    error('showtable:argsize','rowgap option is the wrong size');
end;
if (numel(opt.colgap) == 1),
    colgap = opt.colgap * ones(1,nc);
    colgap([1 nc+1]) = 0;
elseif (length(opt.colgap) == nc-1),
    colgap = [0 opt.colgap(:)' 0];
elseif (length(opt.colgap) == nc),
    colgap = [opt.colgap(:)' 0];
else
    error('showtable:argsize','colgap option is the wrong size');
end;

%get the size of the text and/or information about the plot
set(ax,'Units','points');
h = -1*ones(nr,nc);
width = zeros(nr,nc);
height = zeros(nr,nc);
xl = zeros(nr,nc,2);
yl = zeros(nr,nc,2);
for i = 1:nr,
    for j = 1:nc,
        if (ischar(fmt{i,j}) && ~isempty(fmt{i,j})),
            str = sprintf(fmt{i,j}, X{i,j});
            h(i,j) = text(0,0, str, ...
                'Parent',ax, 'Units','points',...
                textprops{i,j}{:});
            ext1 = get(h(i,j),'Extent');
            width(i,j) = ext1(3) + (colgap(j) + colgap(j+1))/2;
            height(i,j) = ext1(4) + (rowgap(i) + rowgap(i+1))/2;
            xlim{i,j} = [];
            ylim{i,j} = [];
        elseif (isa(fmt{i,j},'function_handle')),            
            h(i,j) = axes('Units','points','Position',[0 0 10 10]);
            feval(fmt{i,j},h(i,j),X{i,j}{:});
            axis(h(i,j),'tight','off');
            xl(i,j,:) = get(h(i,j),'XLim');
            yl(i,j,:) = get(h(i,j),'YLim');
        end;
    end;
end;

pos = get(ax,'Position');

axleft = pos(1);
axbottom = pos(2);
axheight = pos(4);
axwidth = pos(3);

set(ax,'XLim',[0 axwidth], 'YLim',[0 axheight]);

%resize the columns appropriately
if (~isempty(opt.colwidth) && ~strcmp(opt.colwidthmode,'span')),
    opt.colwidthmode = 'manual';
end;
switch lower(opt.colwidthmode),
    case 'even',
        colwidth = max(width(:));
        colwidth = repmat(colwidth,[1 nc]);
    case 'span',
        if (isempty(opt.colwidth)),
            colwidth = (axwidth-sum(colgap)) / nc * ones(1,nc);
        else
            colwidth = (axwidth-sum(colgap)) / nc * opt.colwidth / sum(opt.colwidth);
        end;
    case 'tight',
        colwidth = max(width,[],1);
        if (iscolnames)
            zerowidth = all(width(2:end,:) == 0, 1);
        else
            zerowidth = all(width == 0, 1);
        end;
        if (any(zerowidth)),
            remainingwidth = axwidth - sum(colgap) - sum(colwidth(~zerowidth));
            colwidth(zerowidth) = remainingwidth / sum(zerowidth);
        end;
    case 'manual',
        colwidth = opt.colwidth;
        if (length(colwidth) == 1)
            colwidth = colwidth * ones(1,nc);
        end;
    otherwise,
        error('Unrecognized column width option');
end;

totalwidth = sum(colwidth) + sum(colgap);
if (totalwidth > axwidth),
    colwidth = colwidth * axwidth/totalwidth;
    totalwidth = sum(colwidth);
end;

%and the rows
if (~isempty(opt.rowheight) && ~strcmp(opt.rowheightmode,'span')),
    opt.rowheightmode = 'manual';
end;
rowheight = zeros(nr,1);
if (isempty(opt.colnamesheight))
    r1 = 1;
    nr1 = nr;
else
    rowheight(1) = opt.colnamesheight;
    r1 = 2;
    nr1 = nr - 1;
end;
switch lower(opt.rowheightmode),
    case 'even',
        rowheight = max(flatten(height(r1:end,:)));
        rowheight = repmat(rowheight,[nr 1]);
    case 'span',
        if (isempty(opt.rowheight)),
            rowheight(r1:end) = (axheight - sum(rowgap) - rowheight(1)) / nr1 * ones(nr1,1);
        else
            rowheight(r1:end) = (axheight - sum(rowgap) - rowheight(1)) / nr1 * ...
                rowheight(r1:end) / sum(rowheight(r1:end));
        end;
    case 'tight',
        rowheight(r1:end) = max(height(r1:end,:),[],2);
        if (isrownames),
            zeroheight = all(height(:,2:end) == 0, 2);
        else
            zeroheight = all(height == 0, 2);
        end;
        if (any(zeroheight)),
            remainingheight = axheight - sum(rowgap) - sum(rowheight);
            rowheight(zeroheight) = remainingheight / sum(zeroheight);
        end;
    case 'manual',
        rowheight = opt.rowheight;
        if (length(rowheight) == 1)
            rowheight = rowheight * ones(nr,1);
        end;
    otherwise,
        error('Unrecognized row height option');
end;

totalheight = sum(rowheight);
if (totalheight > axheight),
    ch = cumsum(rowheight);
    showrows = last(ch < axheight);
    totalheight = ch(showrows);
else
    showrows = nr;
end;

%now sort out the x and y limits in the plots
y1 = axheight - cumsum(rowheight(showrows:-1:1)) - cumsum(rowgap(showrows:-1:1));
x1 = [0 cumsum(colwidth(1:end-1))] + cumsum(colgap(1:nc));

[x,y] = meshgrid(x1,y1);

islinkxrows = cellfun(@(x) (ischar(x) && strcmpi(x,'linkrow')), xlim);
islinkxcols = cellfun(@(x) (ischar(x) && strcmpi(x,'linkcol')), xlim);
islinkyrows = cellfun(@(y) (ischar(y) && strcmpi(y,'linkrow')), ylim);
islinkycols = cellfun(@(y) (ischar(y) && strcmpi(y,'linkcol')), ylim);
    
%x and y limits
for i = 1:nr,
    for j = 1:nc,
        if (islinkxrows(i,j)),
            xl1 = [min(xl(i,islinkxrows(i,:),1)) max(xl(i,islinkxrows(i,:),2))];
            xl(i,j,:) = xl1;
        elseif (islinkxcols(i,j)),
            xl1 = [min(xl(islinkxcols(:,j),j,1)) max(xl(islinkxcols(:,j),j,2))];
            xl(i,j,:) = xl1;
        elseif (~isempty(xlim{i,j}) && isnumeric(xlim{i,j}) && (length(xlim{i,j}) == 2)),
            xl(i,j,:) = xlim{i,j};
        end;
        
        if (islinkyrows(i,j)),
            yl1 = [min(yl(i,islinkyrows(i,:),1)) max(yl(i,islinkyrows(i,:),2))];
            yl(i,j,:) = yl1;
        elseif (islinkxcols(i,j)),
            yl1 = [min(yl(islinkycols(:,j),j,1)) max(yl(islinkycols(:,j),j,2))];
            yl(i,j,:) = yl1;
        elseif (~isempty(ylim{i,j}) && isnumeric(ylim{i,j}) && (length(ylim{i,j}) == 2)),
            yl(i,j,:) = ylim{i,j};
        end;
    end;
end;

%and move everything into the correct position
for i = 1:showrows,
    for j = 1:nc,
        if (~ishandle(h(i,j))),
            continue;
        end;
        switch (get(h(i,j),'Type')),
            case 'text',
                set(h(i,j),'Position',[x(i,j) + offx(i,j)*(colwidth(j)-(colgap(j)-colgap(j+1))/2) ...
                    y(i,j) + offy(i,j)*(rowheight(i)-(rowgap(i)-rowgap(i+1))/2)]);
            case 'axes',
                h1 = findobj(h(i,j),'-property','XData');

                for k = 1:length(h1),
                    xd = get(h1(k),'XData');
                    yd = get(h1(k),'YData');
                
                    scalex = (colwidth(j)-(colgap(j)-colgap(j+1))/2)/diff(xl(i,j,:));
                    scaley = (rowheight(i)-(rowgap(i)-rowgap(i+1))/2)/diff(yl(i,j,:));
                    xd = (xd-xl(i,j,1))*scalex + x(i,j);
                    yd = (yd-yl(i,j,1))*scaley + y(i,j);
                    
                    set(h1(k),'Parent',ax, 'XData',xd, 'YData',yd);                            
                end;
                
                h1 = findobj(h(i,j),'-property','BaseValue');
                for k = 1:length(h1),
                    bv = get(h1(k),'BaseValue');
                    bv = (bv-yl(i,j,1))*scaley + y(i,j);
                    set(h1(k),'BaseValue',bv);
                end;
                
                delete(h(i,j));
                h(i,j) = 0;
        end;
    end; 
end;
set(h(showrows+1:end,:),'Visible','off');

set(ax,'Visible','off');
axes(origax);

%------------------------------------------------------------
function [ok,props] = parseheaderfont(len, props, basefont)

if (isempty(props)),
    props = repmat(basefont,[1 len]);
elseif (ischar(props)),
    props = cell(1,len);
    [props{:}] = {'FontName',props};
elseif (iscellstr(props)),
    if ((length(props) == len) && ...
            all(~ismember(props,{'FontName','FontAngle','FontWeight','FontSize'}))),
        props = cellfun(@(x) ({'FontName',x}), props, 'UniformOutput',false);
    else
        props = repmat({props},[1 len]);
    end;
elseif (iscell(props) && (length(props) == len)),
    % everything's OK
else
    ok = false;
    return;
end;

for i = 1:length(props),
    %look for font properties specified in the base font, but not in the
    %current one
    
    %this looks for property names
    isbase = ~ismember(basefont(1:2:end),props{i}(1:2:end));
    
    %now double the length of isbase by replicated each element, so that we
    %get the property values corresponding to each name
    isbase = [isbase; isbase];
    isbase = isbase(:);
    
    props{i} = [basefont(isbase) props{i}];
end;

ok = true;



