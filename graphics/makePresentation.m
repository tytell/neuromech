function makePresentation(varargin)
%MAKEPRESENTATION    Rescales a figure for a slide presentation.
% makePresentation(opts...)
%    or makePresentation(fig,opts...)
% Attempts to rescale figures so that the output is appropriate for
% a presentation.  If no options are given, then it will display a dialog
% asking for the assorted options.
% 
% Adds a button to the rescaled figure to export the figure.  Uses
% EXPORT_FIG from Oliver Woodford.
%
% Standard options:
%   'font', font - Name of the font to use (default: 'Times New Roman')
%   'secondaryfont', font - Name of the secondary font to use for small
%     text, like numbers on axes.  If empty, uses the same as the primary
%     font.
%   'fontsizes', [big small] - Two element vector with the point sizes for
%     small text (numbers on axes) and big text (axes labels)
%   'linescale', m - Multiply all line widths by m.  (default: 2)
%   'minlinewidth', m - Minimum line width in points, after scaling.  (default: 1pt)
%   'markerscale', m - Multiply marker sizes by m.  (default: 2)
%   'minmarkersize', m - Minimum marker size in points, after scaling.
%     (default: 12pt)
%   'bgcolor', col - Usually 'white' or 'black'.  But could be any color
%     string that Matlab supports.  (default: 'white')
%   'replaceblack' - Change black colors to white.  Useful for a black
%     background.
%
% Advanced options:
%   'quiet' - Don't display the dialog.
%   'copyfig', b - Copy the figure before rescaling.  (default: true)
%   'width', w - Width of the output figure in inches.  (default: 10.5in)
%
% See also MAKEPUBLICATION and EXPORT_FIG

% Mercurial revision hash: $Revision$ $Date$
% See https://bitbucket.org/tytell/matlab_graphics/overview
% Copyright (c) 2011, Eric Tytell <tytell at jhu dot edu>

opt.font = 'Times New Roman';
opt.secondaryfont = '';
opt.fontsizes = [24 16];
opt.minlinewidth = 1;
opt.linescale = 2;
opt.axlinewidth = 2;
opt.minmarkersize = 12;
opt.markerscale = 2;
opt.quiet = false;
opt.copyfig = true;
opt.width = 10.5;
opt.bgcolor = 'white';
opt.replaceblack = false;
opt.rescale = true;

if ((nargin >= 1) && (numel(varargin{1}) == 1) && ishandle(varargin{1}))
    fig = varargin{1};
    i = 2;
else
    fig = gcf;
    i = 1;
end;

opt = parsevarargin(opt,varargin(i:end),i);

if (~opt.quiet)
    promptdefault = {'Font name:', opt.font; ...
        'Font size:', num2str(opt.fontsizes(1)); ...
        'Secondary font name:', opt.secondaryfont; ...
        'Secondary font size:', num2str(opt.fontsizes(2)); ...
        'Minimum line thickness (pt):', num2str(opt.minlinewidth);
        'Multiply line thickness:', num2str(opt.linescale);
        'Minimum marker size (pt):', num2str(opt.minmarkersize);
        'Multiply marker size:', num2str(opt.markerscale);
        'Background color:', opt.bgcolor};
    
    figdefs = inputdlg(promptdefault(:,1), 'Figure defaults', 1, promptdefault(:,2));

    opt.font = figdefs{1};
    opt.fontsizes(1) = str2double(figdefs{2});
    opt.secondaryfont = figdefs{3};
    opt.fontsizes(2) = str2double(figdefs{4});
    opt.minlinewidth = str2double(figdefs{5});
    opt.linescale = str2double(figdefs{6});
    opt.minmarkersize = str2double(figdefs{7});
    opt.markerscale = str2double(figdefs{8});
    opt.bgcolor = figdefs{9};
    
    if (strcmpi(opt.bgcolor,'black') || ...
            strcmpi(opt.bgcolor,'k'))
        opt.replaceblack = true;
    end;
end;

if (isempty(opt.secondaryfont))
    opt.secondaryfont = opt.font;
end;

if (opt.copyfig)
    if (strcmp(get(fig,'Tag'),'presentationfigure')),
        figcopy = fig;
        fig = get(figcopy,'UserData');
        delete(figcopy);
    end;
    
    figcopy = findobj('Tag','presentationfigure','UserData',fig);
    if (~isempty(figcopy))
        delete(figcopy);
    end;
    
    figcopy = figure;
    pos = get(fig,'Position');
    units = get(fig,'Units');
    set(figcopy,'Units',units,'Position',pos);
    
    children = get(fig,'Children');
    for i = length(children):-1:1,
        copyobj(children(i),figcopy);
    end;
    
    cmap = get(fig,'Colormap');
    set(figcopy,'Colormap',cmap);
    amap = get(fig,'Alphamap');
    set(figcopy,'Alphamap',amap);
    
    %figcopy = copyobj(fig,0);
    set(figcopy,'Name','presentation copy','UserData',fig,'Tag','presentationfigure', ...
        'Color',opt.bgcolor);
    fig = figcopy;
    
    exportbutt = uicontrol(fig, 'Style','pushbutton','String','Export', ...
        'Position',[20 20 100 35], 'Callback',{@exportpresfig,fig});
end;

allaxes = findobj(fig,'Type','axes');
if (opt.rescale)
    set(fig,'Units','inches');
    
    %set all axes on the figure to have normalized units, so that they scale
    %with the figure when we resize it
    set(allaxes,'Units','normalized','ActivePositionProperty','OuterPosition', ...
        'Clipping','off');
    
    %and resize the figure
    %now set the figure position
    figpos = get(fig,'Position');
    
    scale = opt.width/figpos(3);
    figpos(3) = opt.width;
    
    h0 = figpos(4);
    figpos(4) = h0*scale;
    figpos(2) = figpos(2)+h0 - figpos(4);
    
    set(fig,'Position',figpos);
    
    drawnow;
    %for some reason, a brief pause here seems to help Matlab get the figure
    %sizes right...
    pause(0.05);
end;

%check for markers that are too small
markers = findobj(fig,'Type','line','-not','Marker','none');
if (~isempty(markers)),
    markersz = get(markers,'MarkerSize');
    if (iscell(markersz))
        markersz = cat(1,markersz{:});
    end;

    markersz = markersz * opt.markerscale;
    for i = 1:length(markersz),
        if (markersz(i) < opt.minmarkersize),
            set(markers(i),'MarkerSize',opt.minmarkersize);
        else
            set(markers(i),'MarkerSize',markersz(i));
        end;
    end;
end;

%check for lines that are too thin
%first on the edges of markers
markeredges = findobj(fig,'Type','line','LineStyle','none',...
    '-not','MarkerEdgeColor','none');
markeredgewidth = get(markeredges,'LineWidth');
if (iscell(markeredgewidth))
    markeredgewidth = cat(1,markeredgewidth{:});
end;

set(markeredges(markeredgewidth < opt.minlinewidth),...
    'LineWidth',opt.minlinewidth);

%then on the lines themselves
lines = findobj(fig,'Type','line','-not','LineStyle','none');
if (~isempty(lines)),
    lnwidth = get(lines,'LineWidth');
    if (iscell(lnwidth))
        lnwidth = cat(1,lnwidth{:});
    end;

    lnwidth = lnwidth * opt.linescale;
    for i = 1:length(lnwidth),
        if (lnwidth(i) < opt.minlinewidth),
            set(lines(i),'LineWidth',opt.minlinewidth);
        else
            set(lines(i),'LineWidth',lnwidth(i));
        end;
    end;
end;

%also on axes
axlnwidth = get(allaxes,'LineWidth');
if (iscell(axlnwidth))
    axlnwidth = cat(1,axlnwidth{:});
end;

axlnwidth = opt.axlinewidth;
for i = 1:length(axlnwidth),
    set(allaxes(i),'LineWidth',axlnwidth(i));
end;

%now get the axes labels
labels = zeros(length(allaxes),3);
for i = 1:length(allaxes),
    labels(i,1) = get(allaxes(i),'XLabel');
    labels(i,2) = get(allaxes(i),'YLabel');
    labels(i,3) = get(allaxes(i),'Title');
end;

%set the tick label font and size
set(allaxes,'FontName',opt.font, 'FontSize',opt.fontsizes(2));

good = ishandle(labels);
set(labels(good),'FontName',opt.font, 'FontSize',opt.fontsizes(1));

othertext = findobj(fig,'Type','text');
othertext = othertext(~ismember(othertext,labels(:)));
set(othertext,'FontName',opt.font,'FontSize',opt.fontsizes(1),'FontWeight','bold');

%change black into white if necessary
if (opt.replaceblack)
    colopts = {'Color','MarkerFaceColor','MarkerEdgeColor','EdgeColor','FaceColor', ...
        'XColor','YColor','ZColor'};
    for col = colopts,
        col = col{:};
        black = findobj(get(fig,'Children'),col,'k', '-or', col,'black', '-or', col,[0 0 0]);
        set(black,col,'w');
    end;
    
    black = findobj(labels(:),'Color','k', '-or','Color','black', '-or','Color',[0 0 0]);
    set(black,'Color','w');
    
    whiteaxes = findobj(allaxes,'Color','w', '-or','Color','white', '-or','Color',[1 1 1], ...
        '-depth',0);
    set(whiteaxes,'Color','k');
end;

%make sure that Matlab won't rescale things when we print
set(fig, 'PaperPositionMode','manual');
if (opt.rescale)
    set(fig,'PaperPosition',[0.5 10.5-figpos(4) figpos(3) figpos(4)]);
end;
set(allaxes, 'XTickMode','manual','YTickMode','manual','ZTickMode','manual');




function exportpresfig(hobj,event, fig)     %#ok

set(hobj, 'Visible','off');
[fn,pn] = uiputfile({'*.pdf';'*.jpg';'*.png';'*.tif';'*.eps';'*.*'}, 'Export figure');

export_fig(fullfile(pn,fn),fig);

set(hobj, 'Visible','on');
