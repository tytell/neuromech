function saveallfigs(varargin)
% function saveallfigs()
%  or      saveallfigs(figs)
%  or      saveallfigs(basename)
%  or      saveallfigs(basename,figs)
%
% Saves all figures, or a selection of figures.  Uses basename as the base
% name for the figures (default 'fig'), then appends a number corresponding
% to the figure number.

% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2011, Eric Tytell <tytell at jhu dot edu>

isbasename = false;
switch nargin,
    case 0,
        figs = findobj('Type','figure');
        basename = 'fig';
        
    case 1,
        if (isnumeric(varargin{1})),
            figs = varargin{1};
            basename = 'fig';
        else
            basename = varargin{1};
            figs = findobj('Type','figure');
            isbasename = true;
        end;
        
    case 2,
        figs = varargin{1};
        basename = varargin{2};
        isbasename = true;
end;

for i = 1:length(figs),
    if (~isbasename),
        name = get(figs(i),'Name');
    end;
    if (isempty(name)),
        name = strcat(basename,num2str(figs(i)));
    end;
        
    fprintf('%s...',name);
    saveas(figs(i),name,'fig');
end;
    
   