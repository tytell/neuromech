function makePublication(varargin)
%MAKEPUBLICATION  Rescales a figure for publication in particular journals.
% makePublication(fig, opts...)
% Attempts to rescale a figure for publication in particular journals.  Can
% set line widths, fonts, figure width, and margin.  Tries to duplicate the
% figure so that the original is not affected, but sometimes this causes
% problems.
% 
% Options:
%   'copyfig', true/false - Copy the figure
%   'ticklen', points - Set an absolute tick length, rather than scaling by axes
%      size.

% Mercurial revision hash: $Revision$ $Date$
% See https://bitbucket.org/tytell/matlab_graphics/overview
% Copyright (c) 2011, Eric Tytell <tytell at jhu dot edu>

opt.copyfig = true;
opt.ticklen = 0.020;        % inches

if ((nargin >= 1) && (numel(varargin{1}) == 1) && ...
        ishandle(varargin{1}))
    fig = varargin{1};
    p = 2;
else
    fig = gcf;
    p = 1;
end;
    
opt = parsevarargin(opt,varargin(p:end),p-1);

journaldefaults = {
    'JEB 1 column', {'Times New Roman','10','','9','3.48','0.02'}; ...
    'JEB 2 column', {'Times New Roman','10','','10','7.2','0.02'}; ...
    'J Neurophys narrow', {'Arial','10','Times New Roman','9','3.5','0.02'}; ...
    'J Neurophys medium', {'Arial','10','Times New Roman','9','5','0.02'}; ...
    'J Neurophys wide', {'Arial','10','Times New Roman','9','7','0.02'}; ...
    'J Neurosci 1 col', {'Arial','10','','9','3.345','0.02'}; ...
    'J Neurosci 1.5 col', {'Arial','10','','9','4.56','0.02'}; ...
    'J Neurosci 2 col', {'Arial','10','','9','6.93','0.02'}; ...
    'Exp Fluids narrow', {'Times New Roman','10','','9','3.38','0.02'}; ...
    'PLOS 1 column', {'Arial','10','','9','3.32','0.02'}; ...
    'PLOS 1.5 column', {'Arial','10','','9','4.92','0.02'}; ...
    'PLOS 2 column', {'Arial','10','','9','6.89','0.02'}; ...
    'ICB 1 column', {'Arial','10','','9','3.46','0.02'}; ...
    'ICB 2 column', {'Arial','10','','9','7.09','0.02'}; ...
    'PNAS 1 column', {'Arial','10','','9','3.42','0.02'}; ...
    'PNAS 1.5 column', {'Arial','10','','9','4.5','0.02'}; ...
    'PNAS 2 column', {'Arial','10','','9','7','0.02'}; ...
    'Elsevier 0.5 col', {'Helvetica','10','','9','1.535','0.02'};
    'Elsevier 1 column', {'Helvetica','10','','9','3.307','0.02'};
    'Elsevier 1.5 column', {'Helvetica','10','','9','5.079','0.02'};
    'Elsevier 2 column', {'Helvetica','10','','9','6.85','0.02'};
    'Poster large', {'Helvetica','36','','21','11','0'}; ...
    'Custom',{'','','','','',''}};

[sel,ok] = listdlg('PromptString','Select journal style',...
    'SelectionMode','single','ListString',journaldefaults(:,1));

if (~ok)
    return;
end;

promptdefault = {'Primary font name:', ''; ...
                 'Primary font size:', ''; ...
                 'Secondary font name:', ''; ...
                 'Secondary font size:', ''; ...
                 'Figure width (in):', ''; ...
                 'Margin (in):', ''};
[promptdefault{:,2}] = deal(journaldefaults{sel,2}{:});

figdefs = inputdlg(promptdefault(:,1), 'Figure defaults', 1, promptdefault(:,2));

fontname1 = figdefs{1};
fontsize1 = str2double(figdefs{2});
fontname2 = figdefs{3};
if (isempty(fontname2))
    fontname2 = fontname1;
end;
fontsize2 = str2double(figdefs{4});
figwidth = str2double(figdefs{5});
figmargin = str2double(figdefs{6});

if (opt.copyfig)
    if (strcmp(get(fig,'Tag'),'publicationfigure')),
        figcopy = fig;
        fig = get(figcopy,'UserData');
        delete(figcopy);
    end;
    
    figcopy = findobj('Tag','publicationfigure','UserData',fig);
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
    set(figcopy,'Name','publication copy','UserData',fig,'Tag','publicationfigure');
    fig = figcopy;

    uicontrol(fig, 'Style','pushbutton','String','Export', ...
        'Position',[20 20 100 35], 'Callback',{@exportpubfig,fig});
end;

set(fig,'Units','inches');

%and resize the figure
%now set the figure position
figpos = get(fig,'Position');

scale = figwidth/figpos(3);
figpos(3) = figwidth;

h0 = figpos(4);
figpos(4) = h0*scale;
figpos(2) = figpos(2)+h0 - figpos(4);

set(fig,'Position',figpos);

%adjust the text sizes
allaxes = findobj(fig,'Type','axes');

labels = zeros(length(allaxes),3);
for i = 1:length(allaxes),
    labels(i,1) = get(allaxes(i),'XLabel');
    labels(i,2) = get(allaxes(i),'YLabel');
    labels(i,3) = get(allaxes(i),'Title');
end;

%set the tick label font and size
set(allaxes,'FontName',fontname2, 'FontSize',fontsize2, 'LineWidth',1);

good = ishandle(labels);
set(labels(good),'FontName',fontname1, 'FontSize',fontsize1);

othertext = findobj(fig,'Type','text');
othertext = othertext(~ismember(othertext,labels(:)));
set(othertext,'FontName',fontname1,'FontSize',fontsize1);

set(allaxes,'LineWidth',1);
set(allaxes,'Units','normalized','ActivePositionProperty','OuterPosition');

drawnow;
pause(0.05);

%the axes should be about the right size
%get the boundaries of all the axes
bound = [Inf Inf 0 0];
for i = 1:length(allaxes)
    set(allaxes(i),'Units','inches');
    
    tight = get(allaxes(i),'TightInset');
    axpos = get(allaxes(i),'Position');
    
    bound(1) = min(bound(1),axpos(1)-tight(1));
    bound(2) = min(bound(2),axpos(2)-tight(2));
    %bound 3 and 4 are right and top, not width and height
    bound(3) = max(bound(3),axpos(1)+axpos(3)+tight(3));
    bound(4) = max(bound(4),axpos(2)+axpos(4)+tight(4));
    
    ticklen1 = get(allaxes(i),'TickLength');
    ticklen1(1) = opt.ticklen / max(axpos(3:4));
    set(allaxes(i),'TickLength',ticklen1);
    
    htxt = findobj(allaxes(i),'Type','text');
    for j = 1:length(htxt),
        units = get(htxt(j),'Units');
        set(htxt(j),'Units','inches');
        txtpos = get(htxt(j),'Extent');
        txtpos([1 2]) = txtpos([1 2]) + axpos([1 2]);
        set(htxt(j),'Units',units);
        
        bound(1) = min(bound(1),txtpos(1));
        bound(2) = min(bound(2),txtpos(2));
        %bound 3 and 4 are right and top, not width and height
        bound(3) = max(bound(3),txtpos(1)+txtpos(3));
        bound(4) = max(bound(4),txtpos(2)+txtpos(4));
    end;
end;
%turn bound(3) and (4) back to width and height
width = bound(3) - bound(1);
height = bound(4) - bound(2);

scale = (figwidth-2*figmargin) / width;

%and scale the axes
for i = 1:length(allaxes),
    pos = get(allaxes(i), 'Position');
    pos(3) = pos(3)*scale;
    pos(4) = pos(4)*scale;
    pos(1) = (pos(1) - bound(1))*scale + figmargin;
    pos(2) = height*scale - figmargin - (bound(4) - pos(2))*scale;
    
    set(allaxes(i), 'Position',pos);
end;

set(fig, 'PaperPositionMode','manual',...
    'PaperPosition',[0.5 10.5-figpos(4) figpos(3) figpos(4)], 'Color','w');
set(allaxes, 'XTickMode','manual','YTickMode','manual','ZTickMode','manual');



function exportpubfig(hobj,event, fig)     %#ok

set(hobj, 'Visible','off');
[fn,pn] = uiputfile({'*.pdf';'*.jpg';'*.png';'*.tif';'*.eps';'*.*'}, 'Export figure');

export_fig(fullfile(pn,fn),fig);

set(hobj, 'Visible','on');
