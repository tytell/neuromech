function h = mplot(varargin)
% h = mplot(x,y,multilinespec)
%
% Works more or less like plot, but it can take multiple line spec
% parameters for multiple plots.  Separate different linespecs by '|' or in
% cell arrays.  Adds the option 'f' to the standard linespec parametrs, which means
% that the marker will be filled.  Standard options for plot can be passed with multiple
% arguments.
%
% Example:
%   plot(x,y,'ks|ro')
%
%   Plots x and y, alternating between black squares and red circles as the
%   markers.
%
%   plot(x,y,'of','Color',[0.7 0.7 0.7; 0.5 0 0])
%
%   Plots x and y, alternating between gray circles and dark red circles,
%   both filled.
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell

if ((numel(varargin{1}) == 1) && ishandle(varargin{1}) && ...
    strcmp(get(varargin{1},'Type'),'axes')),
    ax = varargin{1};
    p = 2;
else
    ax = gca;
    p = 1;
end;

%matches a linespec or multilinespec (will also match improperly formed linespecs...)
lspecchars = '-:.+o*.xsd^v<>phrgbcmykwf|';

ischararg = cellfun(@ischar,varargin);
islinespec = false(size(varargin));
islinespec(ischararg) = cellfun(@(x) (all(ismember(x,lspecchars))), varargin(ischararg));

if (any(cellfun(@ndims, varargin(~ischararg)) > 2)),
    error('mplot:toomanydims','Cannot handle 3D matrices');
end;

xyl = {};

%process the x,y,linespec options
a = 1;
done = false;
while (~done && (p <= nargin)),
    if ((p+1 <= nargin) && isnumeric(varargin{p}) && isnumeric(varargin{p+1})),
        x1 = varargin{p};
        y1 = varargin{p+1};
        
        if (any(size(x1) == 1)),
            %check for vector x and matrix y conditions and make sure the
            %dimension with the same length is 1
            if (size(y1,2) == numel(x1)),
                y1 = y1';
            end;
            x1 = x1(:);
            x1 = x1(:,ones(1,size(y1,2)));
        elseif (any(size(y1) == 1)),
            %check for matrix x and vector y
            if (size(x1,2) == numel(y1)),
                x1 = x1';
            end;
            y1 = y1(:);
            y1 = y1(:,ones(1,size(x1,2)));
        end;
        
        n = size(y1,2);
        
        xyl1 = cell(3,n);
        
        %separate out the columns of x and y into cells
        xyl1(1,:) = mat2cell(x1, size(x1,1), ones(1,n));
        xyl1(2,:) = mat2cell(y1, size(x1,1), ones(1,n))';
    
        p = p+2;
        
        %look for a linespec
        if ((p <= nargin) && islinespec(p)),
            if (iscellstr(varargin{p})),
                %cellstrs are easy
                ls = varargin{p};
            else
                %split a string on |
                ls = regexp(varargin{p},'\|','split');
            end;
            
            %repeat it the appropriate number of times
            nrep = ceil(n/length(ls));
            ls = repmat(ls(:)',[1 nrep]);
            ls = ls(1:n);
            
            %and make that the third row
            xyl1(3,:) = ls;
            
            p = p+1;
        end;
        
        %collect all of the x,y,linespec triples
        xyl{a} = xyl1;
        a = a+1;
    else
        %finish when we don't match an x,y combination
        done = true;
    end;
end;
        
%squish all the x,y,linespec triples together
xyl = cat(2,xyl{:});
n = size(xyl,2);

%look for options
opt = cell(0,n);

while (p <= nargin),   
    if (~ischar(varargin{p})),
        error('mplot:badoption','Unrecognized option at parameter %d',p);
    else
        %options apply to all x,y,linespec parameters
        opt1 = cell(2,n);
        [opt1{1,:}] = deal(varargin{p});
        if (p+1 > nargin),
            error('mplot:missingoptionarg','Option %s requires an argument',varargin{p});
        else
            optval = varargin{p+1};
            
            %check for options with multiple values
            if (iscell(optval)),
                %cell arrays just get repeated to match the number of lines
                optval = repmat(optval(:)',[1 ceil(n/numel(optval))]);
                opt1(2,:) = optval(1:n);
            elseif (ischar(optval) && any(optval == '|')),
                %split character arrays on |
                optval = regexp(optval,'\|','split');
                optval = repmat(optval(:)',[1 ceil(n/numel(optval))]);
                opt1(2,:) = optval(1:n);
            elseif (strcmp(varargin{p},'Color') && isnumeric(optval) && ...
                    any(size(optval) == 3)),
                %deal with the 'Color' option differently, because it has
                %sets of repeated triples
                if (size(optval,2) == 3),
                    optval = optval';
                end;
                optval = repmat(optval,[1 ceil(n/size(optval,2))]);
                optval = optval(:,1:n);
                optval = mat2cell(optval,3,ones(1,n));
                opt1(2,:) = optval;
            elseif (ischar(optval) && ~isempty(strmatch(optval,{'on','off','none','auto'}))),
                %handle single keywords
                [opt1{2,:}] = deal(optval);
            else
                %otherwise, split on single valuse and replicate to match
                %the number of lines
                optval = repmat(optval(:)',[1 ceil(n/size(optval,2))]);
                optval = optval(1:n);
                optval = num2cell(optval);
                opt1(2,:) = optval;
            end;
            p = p+2;
        end;
    
        opt(end+1:end+2,:) = opt1;
    end;
end;

%handle filled symbols
%look for the 'f' option
isfill = cellfun(@(x) (any(x == 'f')), xyl(3,:));
%remove it, because plot doesn't understand it
xyl(3,:) = cellfun(@(x) (x(x ~= 'f')), xyl(3,:), 'UniformOutput',false);

%make sure empty linespecs are of type char
emptylinespec = cellfun(@isempty, xyl(3,:));
[xyl{3,emptylinespec}] = deal('');

%what color are we filling with?
charopt = find(cellfun(@ischar,opt(:,1)));
k = strmatch('Color',opt(charopt,1),'exact');
k = charopt(k);
if (~isempty(k)),
    %use the 'Color' option, if they passed it in
    col = opt(k+1,:);
else
    %otherwise look for the linespec
    col = cellfun(@(x) (x(ismember(x,'ymcrgbwk'))), xyl(3,:),'UniformOutput',false);
    if (all(cellfun(@isempty,col))),
        col0 = get(ax,'ColorOrder');
        ncol = ceil(length(col)/size(col0,1));
        col0 = repmat(col0,[ncol 1]);
        col0 = col0(1:length(col),:);
        col = mat2cell(col0,ones(length(col),1),3);
    end;
end;

%get rid of any existing 'MarkerFaceColor' option
k = strmatch('MarkerFaceColor',opt(charopt,1),'exact');
k = charopt(k);
if (~isempty(k)),
    mfcol = opt(k+1,:);
    opt = opt([1:k-1 k+2:end],:);
    ismfcol = true;
else
    mfcol = cell(1,n);
    ismfcol = false;
end;

[mfcol{cellfun(@isempty,mfcol)}] = deal('auto');

if (any(isfill) && ismfcol),
    warning('mplot:optoverride','Fill option overrides MarkerFaceColor');
end;

%and set the 'MarkerFaceColor' appropriately for the fill
nopt = size(opt,1);
[opt{nopt+1,:}] = deal('MarkerFaceColor');
if (any(isfill)),
    opt(nopt+2,isfill) = col(isfill);
end;

if (any(~isfill)),
    opt(nopt+2,~isfill) = mfcol(~isfill);
end;

%do the plot
h = plot(ax,xyl{:});
for i = 1:length(h),
    %set the options
    set(h(i),opt{:,i});
end;



                    