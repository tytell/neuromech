function figblackbg(fig)
% Sets a figure to have a black background with white axes, for black
% background PowerPoint presentations
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell <tytell at jhu dot edu>

if (nargin == 0),
    fig = gcf;
end;

set(fig,'Color','k');

%find the axes and set them to be transparent with white lines
ax = findobj(fig,'Type','axes');
set(ax,'Color','none');
set(ax,'XColor','w','YColor','w');

%find any black line objects and make them white
ln = findobj(ax,'Type','line','Color','k');
set(ln,'Color','w');
mark = findobj(ax,'Type','line','MarkerFaceColor','k');
set(mark,'MarkerFaceColor','w');
mark = findobj(ax,'Type','line','MarkerEdgeColor','k');
set(mark,'MarkerEdgeColor','w');

%same for black text
txt = get(ax,'XLabel');
txt = cat(1,txt,get(ax,'YLabel'));
txt = cat(1,txt,get(ax,'Title'));
if (iscell(txt)),
    txt = cat(1,txt{:});
end;
txt = cat(1,txt,findobj(fig,'Type','text','Color','k'));
set(txt,'Color','w');
