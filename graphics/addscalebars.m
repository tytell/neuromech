function varargout = addscalebars(varargin)
% function hbar = addscalebars(options...)
% 
% Options (in any order) can be
%   'x','y', or 'xy' - X, Y, or both scale bars
%   'xlen','ylen', followed by the length - Defines the length of the scale
%       bars.  If not specified, will use the spacing of the ticks to define the
%       lengths
%   'units', {xunit,yunit} or units - Appends a text label for the units.
%   axes - Axes to add the scale bar to (default is gca)
%   'off' - Removes existing scale bars
%   'hideaxis', true/false - Hides the existing axis after adding the scale
%       bar
%   'position', [x y] - Location of the bar in normalized coordinates.
%       Default is [1 0.1], which means all the way on the right and
%       slightly up from the bottom
%   'textlocation' - Location of the text relative to the bar.  Can have
%       the following values: 'ctr','center','left','l','right','r',
%       'top','t','bottom','b','leftbottom','lb','rightbottom','rb',
%       'lefttop','lt','righttop',rt'.  You can also just write the value,
%       without the text location option
%
% Examples:
%    addscalebars('x','xlen',2,'units','cm','position',[0.8 -0.1]) - Just
%      x scale bar, 2cm long located near the right below the axes.
%    addscalebars('xy','rightbottom','position',[1 0])

% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell <tytell at jhu dot edu>

opt.bartype = 'xy';                         % type of plot
opt.xlen = [];                              % length of x bar
opt.ylen = [];                              % length of y bar
opt.hideaxis = true;                    % hide the existing axis
opt.position = [1 0.1];
opt.textlocation = 'rb';
opt.units = {};

if ((nargin >= 1) && (numel(varargin{1}) == 1) && ishandle(varargin{1}) && ...
        strcmp(get(varargin{1},'Type'), 'axes'))
    ax = varargin{1};
    p = 2;
else
    ax = gca;
    p = 1;
end;

opt = parsevarargin(opt,varargin(p:end), 'typecheck',false,...
    'multival',{'bartype',{'off','update','x','y','xy'}; ...
                'textlocation',{'ctr','center','left','l',...
                               'right','r','top','t','bottom','b',...
                               'leftbottom','lb',...
                               'rightbottom','rb','lefttop','lt',...
                               'righttop','rt'}},...
                'synonyms',{{'textlocation',{'orientation'}},...
                            {'position',{'pos'}}});
                        
%default lengths of the bars are the lengths between one tick on the plot
if (isempty(opt.xlen)),
    xt = get(ax,'XTick');
    if (isempty(xt)),
        set(ax,'XTickMode','auto');
        xt = get(ax,'XTick');
        set(ax,'XTick',[]);
    end;
    opt.xlen = xt(2) - xt(1);
end;

if (isempty(opt.ylen)),
    yt = get(ax,'YTick');
    if (isempty(yt)),
        set(ax,'YTickMode','auto');
        yt = get(ax,'YTick');
        set(ax,'YTick',[]);
    end;
    opt.ylen = yt(2) - yt(1);
end;

axpos = get(ax,'Position');
xlm = get(ax,'XLim');
ylm = get(ax,'YLim');

%find any existing bars
hbar = findobj(ax,'Tag','scalebar');

%check whether they're x, y, or both
existingbartype = '';
if (~isempty(hbar))
    hln = findobj(hbar,'Type','line');
    xd = get(hln,'XData');
    yd = get(hln,'YData');
    if (diff(yd) == 0),
        existingbartype = 'x';
    elseif (diff(xd) == 0),
        existingbartype = 'y';
    else
        existingbartype = 'xy';
    end;
end;

%to update, just delete the existing bar and redo
if (strcmp(opt.bartype,'update')),
    if (isempty(existingbartype)),
        error('No existing scale bar to update');
    end;
    opt.bartype = existingbartype;
end;

if (~isempty(hbar)),
    delete(hbar);
end;

hbar = zeros(3,1);

isheld = ishold(ax);
hold(ax, 'on');

xctr = diff(xlm) * opt.position(1);
yctr = diff(ylm) * opt.position(2);

switch get(ax,'XDir'),
    case 'normal',
        xctr = xctr + xlm(1);
    case 'reverse',
        xctr = xlm(2) - xctr;
end;
switch get(ax,'YDir'),
    case 'normal',
        yctr = yctr + ylm(1);
    case 'reverse',
        yctr = ylm(2) - yctr;
end;

xalign = 'top';
yalign = 'top';
xoff = 0;
yoff = 0;
xd = 1;
yd = 1;

switch lower(opt.textlocation),
    case {'center','ctr'},
        % do nothing
    case {'left','l'},
        yalign = 'bottom';
        xoff = opt.xlen/2;
    case {'right','r'},
        yalign = 'top';
        xoff = -opt.xlen/2;
    case {'top','t'},
        xalign = 'bottom';
        yoff = -opt.ylen/2;
    case {'bottom','b'},
        xalign = 'top';
        yoff = opt.ylen/2;
    case {'leftbottom','lb'},
        yalign = 'bottom';
        xoff = opt.xlen/2;
        yoff = opt.ylen/2;
        xd = -1;
    case {'rightbottom','rb'},
        xoff = -opt.xlen/2;
        yoff = opt.ylen/2;
    case {'lefttop','lt'},
        xalign = 'bottom';
        yalign = 'bottom';
        xoff = opt.xlen/2;
        yoff = -opt.ylen/2;
        xd = -1;
        yd = -1;
    case {'righttop','rt'},
        yalign = 'bottom';
        xoff = -opt.xlen/2;
        yoff = -opt.ylen/2;
        yd = -1;
end;

switch lower(opt.bartype),
    case 'x',
        xx = [-opt.xlen; opt.xlen]/2 + xctr + xoff;
        xtxt = num2str(opt.xlen);
        if (~isempty(opt.units))
            xtxt = strcat(xtxt,opt.units);
        end;
        
        %plot the bar itself
        hbar(1) = plot(ax, xx,[yctr;yctr], 'k-', 'LineWidth',3, 'Clipping','off',...
            'Tag','scalebar');
        %and the text label
        hbar(2) = text(xctr+xoff,yctr, xtxt, 'Parent',ax, ...
            'HorizontalAlignment','center', 'VerticalAlignment',xalign,...
            'Clipping','off','Tag','scalebar');
        
        %hide the corresponding axis if necessary
        if (opt.hideaxis),
            col = get(get(ax,'Parent'),'Color');
            set(ax, 'XTick',[],'XColor',col);
            hlab = get(ax,'XLabel');
            delete(hlab(ishandle(hlab)));
        end;
        
    case 'y',
        yy = [-opt.ylen; opt.ylen]/2 + yctr + yoff;
        ytxt = num2str(opt.ylen);
        if (~isempty(opt.units)),
            ytxt = strcat(ytxt,opt.units);
        end;
        
        hbar(1) = plot(ax, [xctr;xctr],yy, 'k-', 'LineWidth',3, ...
                       'Clipping','off', 'Tag','scalebar');
        hbar(3) = text(xctr,yctr+yoff, ytxt, 'Parent',ax, ...
            'HorizontalAlignment','center', 'VerticalAlignment',yalign, ...
            'Clipping','off','Rotation',90, 'Tag','scalebar');
        if (opt.hideaxis),
            col = get(get(ax,'Parent'),'Color');
            set(ax, 'YTick',[],'YColor',col);
            hlab = get(ax,'YLabel');
            delete(hlab(ishandle(hlab)));
        end;
        
    case 'xy',
        xx = xd*[-opt.xlen opt.xlen opt.xlen]/2 + xctr + xoff;
        yy = yd*[-opt.ylen -opt.ylen opt.ylen]/2 + yctr + yoff;
        xtxt = num2str(opt.xlen);
        ytxt = num2str(opt.ylen);
        if (iscellstr(opt.units) && (length(opt.units) == 2))
            xtxt = strcat(xtxt,opt.units{1});
            ytxt = strcat(ytxt,opt.units{2});
        end;
        
        hbar(1) = plot(ax, xx,yy, 'k-', 'LineWidth',3, 'Clipping','off',...
            'Tag','scalebar');
        hbar(2) = text(xctr+xoff,yy(1), xtxt, 'Parent',ax, ...
            'HorizontalAlignment','center', 'VerticalAlignment',xalign, ...
            'Clipping','off','Tag','scalebar');
        hbar(3) = text(xx(2),yctr+yoff, ytxt, 'Parent',ax, ...
            'HorizontalAlignment','center', 'VerticalAlignment',yalign, ...
            'Clipping','off','Rotation',90, 'Tag','scalebar');
        if (opt.hideaxis),
            col = get(get(ax,'Parent'),'Color');
            set(ax, 'XTick',[],'XColor',col, 'YTick',[],'YColor',col);
            hlab = [get(ax,'XLabel') get(ax,'YLabel')];
            delete(hlab(ishandle(hlab)));
        end;
        
    case 'off',
        if (opt.hideaxis),
            %make the normal axis reappear if necessary
            switch existingbartype,
                case 'x',
                    set(ax,'XColor','k','XTickMode','auto');
                case 'y',
                    set(ax,'YColor','k','YTickMode','auto');
                case 'xy',
                    set(ax,'XColor','k','XTickMode','auto','YColor','k','YTickMode','auto');
            end;
        end;
end;

if (~isheld),
    hold(ax,'off');
end;

if (nargout == 1)
    varargout = {hbar(hbar ~= 0)};
end;









