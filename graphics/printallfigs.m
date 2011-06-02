function printallfigs(varargin)
% function printallfigs(varargin)
% Prints all open figures, using the printpreview dialog for the first one
% and then using the same settings for all the rest.

% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2011, Eric Tytell <tytell at jhu dot edu>

if ((nargin == 0) || ischar(varargin{1})),
    figs = findobj('Type','figure');
else
    figs = varargin{1};
end;

settings = {'PaperUnits',[]; 'PaperOrientation',[]; 'PaperPosition',[]; ...
    'PaperSize',[]; 'PaperType',[]; 'Renderer',[]};
settings = settings';

printpreview(figs(1));
inputyn('Press return when done');

for i = 1:size(settings,2),
    settings{2,i} = get(figs(1),settings{1,i});
end;

if (nargin > 0),
    k1 = find(cellfun(@ischar,varargin));
    k = strncmp('-d',varargin{k1},2);
    k = k1(k);
else
    k = [];
end;
 
basename = 'fig';
if (~isempty(k)),
    if (ismember(lower(varargin{k}(3:end)), ...
            {'eps','epsc','epsc2','jpeg','pdf','tiff'})),
        if (varargin{end}(1) ~= '-'),
            basename = varargin{end};
            varargin = varargin(1:end-1);
        end;
        isname = true;
    end;
else
    isname = false;
end;

datestr = date;
h = zeros(3);
names = cell(length(figs),1);
for i = 1:length(figs),
    names{i} = get(figs(i),'Name');
    
    h(1) = axes('Parent',figs(i),'Position',[0 0 1 1],'Visible','off');
    if (~isempty(names{i})),
        h(2) = text('Parent',h(1),'Units','normalized','Position',[0.1 0.98], ...
            'HorizontalAlignment','left','VerticalAlignment','top',...
            'String',names{i});
    else
        names{i} = sprintf('%s%d',basename,figs(i));
    end;
    
    h(3) = text('Parent',h(1),'Units','normalized','Position',[0.9 0.98], ...
        'HorizontalAlignment','right','VerticalAlignment','top',...
        'String',datestr);
    
    set(figs(i),settings{:});
    
    printfig = sprintf('-f%d',figs(i));
    if (isname),
        opts = [{printfig} varargin names(i)];
    else
        opts = [{printfig} varargin];
    end;
    
    fprintf('%s...\n',names{i});
    
    print(opts{:});
    delete(h(h ~= 0));
end;
    
   