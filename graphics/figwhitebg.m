function figwhitebg(fig)
% Sets a figure to have a white background with black axes
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell <tytell at jhu dot edu>

if (nargin == 0),
    fig = gcf;
end;

set(fig,'Color','w');

%find the axes and set them to be transparent with white lines
ax = findobj(fig,'Type','axes');
set(ax,'Color','none');
set(ax,'XColor','k','YColor','k');

%find any black line objects and make them white
ln = findobj(ax,'Type','line','Color','w');
set(ln,'Color','k');
mark = findobj(ax,'Type','line','MarkerFaceColor','w');
set(mark,'MarkerFaceColor','k');
mark = findobj(ax,'Type','line','MarkerEdgeColor','w');
set(mark,'MarkerEdgeColor','k');

%same for black text
txt = get(ax,'XLabel');
txt = cat(1,txt,get(ax,'YLabel'));
txt = cat(1,txt,get(ax,'Title'));
if (iscell(txt)),
    txt = cat(1,txt{:});
end;
txt = cat(1,txt,findobj(fig,'Type','text','Color','w'));
set(txt,'Color','k');
