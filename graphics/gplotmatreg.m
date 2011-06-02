function gplotmatreg(X,varargin)
% function gplotmatreg(X,[Y],[group],[names],[mksize,txtsize])
%
% X is a n x m matrix, where n is the number of data points and m is the number
% of variables
% Plots all the variables in X against each other
%
% You can optionally pass a grouping variable (group), names of the variables
% (names), and sizes for the plot points (mksize) and text (txtsize) in points.
% 
% names has to be a cell string
% if you want to pass txtsize, you have to give mksize also, before, in
% that order
%
% Also, various text options:
%   loglog - Do log-log plots
%   logx, logy - Make the x or y axis logarithmic
%   fill - Fill the points
%   noaxis - Do not show the axes
%   noticks - Do not show ticks or numbers on the axes
%
% See also:
%   gplotmatrix (Statistics Toolbox)

% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell <tytell at jhu dot edu>

xmargin = 0.1;
ymargin = 0.1;
nbins = 5;
col = 'bgrcmyk';
ncol = 7;
sym = 'o+sd^v><*.';
nsym = 10;
isaxis = 1;
isfill = 0;
isticks = 1;
logx = 0;
logy = 0;

mksize = [];
txtsize = [];
gap = [];
Y = [];
ynames = [];
names = [];
group = [];

N = size(X,1);
n = size(X,2);

for i = 1:length(varargin),
    if ((size(varargin{i},1) == N) & (size(varargin{i},2) > 1)),
        Y = varargin{i};
    elseif (length(varargin{i}) == N),
        group = shiftdim(varargin{i});
    elseif (iscellstr(varargin{i}) & (length(varargin{i}) == n) & ...
            isempty(names)),
        names = varargin{i};
    elseif (iscellstr(varargin{i}) & (length(varargin{i}) == size(Y,2))),
        ynames = varargin{i};
    elseif (isnumeric(varargin{i}) & (length(varargin{i}) == 1)),
        if (isempty(mksize)),
            mksize = varargin{i};
        elseif (isempty(txtsize)),
            txtsize = varargin{i};
        else
            gap = varargin{i};
        end;
    elseif (ischar(varargin{i})),
        switch (lower(varargin{i})),
         case 'fill',
          isfill = 1;
         case 'noaxis',
          isaxis = 0;
         case 'noticks',
          isticks = 0;
         case 'loglog',
          logx = 1;
          logy = 1;
         case 'logx',
          logx = 1;
         case 'logy',
          logy = 1;
         otherwise
          error(sprintf('Unrecognized option %s\n',varargin{i}));
        end;
    else
        error(sprintf('Unknown parameter #%d\n',i+1));
    end;
end;

if (isempty(mksize))
    mksize = 4;
end;
if (isempty(txtsize))
    txtsize = 10;
end;
if (isempty(gap))
    gap = 0;
end;

if (~isempty(group)),
    k = 1;
    gp = double(group);
    for i = 1:length(gp),
        if (isfinite(gp(i))),
            gpk{k} = find(gp == gp(i));
            if (isnumeric(group(i))),
                gpnames{k} = num2str(group(i));
            else
                gpnames{k} = group(i);
            end;
            
            gp(gpk{k}) = NaN;
            k = k+1;
        end;
    end;
    ngp = k-1;
else
    ngp = 1;
    gpk{1} = 1:N;
    gpnames = [];
end;

clf;

[i,j] = find(isnan(X));
XX = X;
if (~isempty(i)),
    meanX = nanmean(X);
    XX(sub2ind(size(XX),i,j)) = meanX(j);
end;
r = corrcoef(XX);

if (isempty(Y))
    nx = n;
    ny = n;
else
    nx = size(X,2);
    ny = size(Y,2);
end;

w = (1 - 2*xmargin)/nx;
h = (1 - 2*ymargin)/ny;

hh = repmat(NaN,[ny nx]);
for i = ny:-1:1,
    for j = 1:nx,
        x = (j-1)*w + xmargin;
        y = 1 - i*h - ymargin;
        
        hh(i,j) = axes('Position',[x y w-gap h-gap],...
                       'FontSize',txtsize);
	
        if ((i == j) & isempty(Y)),
            if (ngp > 1),
                hold on;
                for k = 1:ngp,
                    x = X(:,i);
                    [freq1,xout1] = hist(x(gpk{k}),nbins);
                    plot(xout1,freq1,'-','Color',col(mod(k,ncol)+1));
                end;
                hold off;
            else
                x = X(:,i);
                [freq,xout] = hist(x,nbins);
                bar(xout,freq);
            end;

            dx = 0.05*range(x);
            axis([nanmin(x)-dx nanmax(x)+dx 0 1]);
            axis auto-y;
            
            set(gca,'XTickLabel','','YTickLabel','');
        else
            if (isempty(Y))
                y = X(:,i);
            else
                y = Y(:,i);
            end;
            x = X(:,j);

            xr = x;
            if (logx),
                xr(xr <= 0) = NaN;
                xr = log(x);
            end;
            yr = y;
            if (logy),
                yr(yr <= 0) = NaN;
                yr = log(y);
            end;
                
            [b,bint,r,rint,stats] = regress(yr,[xr ones(size(X,1),1)]);

            hold on;
            for k = 1:ngp,
                if (isfill)
                    mfc = col(mod(k,ncol));
                else
                    mfc = 'none';
                end;
                
                plot(x(gpk{k}),y(gpk{k}),'Color',col(mod(k,ncol)+1),...
                     'Marker',sym(mod(k,nsym)+1),'MarkerSize',mksize,...
                     'LineStyle','none','MarkerFaceColor',mfc);
            end;
            if (sign(bint(1,1)) == sign(bint(1,2))),
                mx = linspace(min(x),max(x),10);
                if (logx),
                    my = polyval(b,log(mx));
                    my = exp(my);
                else
                    my = polyval(b,mx);
                end;
                plot(mx,my,'k-');
            end;
            hold off;
            set(gca,'XTickLabel','','YTickLabel','');
            
            text(0.02,0.98,sprintf('%.3f',stats(1)),'Units','normalized','VerticalAlignment','top');
            
            if (~isticks),
                dx = 0.05*range(x(isfinite(x)));
                dy = 0.05*range(y(isfinite(y)));
                
                axis([nanmin(x)-dx nanmax(x)+dx nanmin(y)-dy nanmax(y)+dy]);
            end;
        end;
        
        if (logx),
            set(gca,'XScale','log');
        end;
        if (logy),
            set(gca,'YScale','log');
        end;
    end;
end;

if (~isempty(gpnames)),
    legend(gpnames{:});
end;

for i = 1:ny,
    if ((i ~= 1) & isticks),
        set(hh(i,1), 'YTickLabelMode','auto');
    end;
    
    nm = [];
    if (isempty(ynames))
        if (~isempty(names)),
            nm = names{i};
        end;
    else
        nm = ynames{i};
    end;
    
    if (~isempty(nm)),
        if (isaxis)
            set(get(hh(i,1),'YLabel'),'String',nm,'FontSize',txtsize);
        else
            set(gcf,'CurrentAxes',hh(i,1));
            text(0,0.5,nm,'FontSize',txtsize,'Rotation',90,'Units','normalized',...
                 'HorizontalAlignment','center','VerticalAlignment','bottom');
        end;
    end;
end;

for j = 1:n,
    if (isticks)
        set(hh(end,j), 'XTickLabelMode','auto');
    end;
    if (~isempty(names)),
        if (isaxis)
            set(get(hh(end,j),'XLabel'),'String',names{j},'FontSize',txtsize);
        else
            set(gcf,'CurrentAxes',hh(end,j));
            text(0.5,0,names{i},'FontSize',txtsize,'Units','normalized',...
                 'HorizontalAlignment','center','VerticalAlignment','top');
        end;
    end;
    
    if (isempty(Y)),
        set(hh(j,j),'XLim', get(hh(end,j),'XLim'));
        set(hh(j,j),'YLimMode','auto');
    end;
end;

if (~isaxis),
    for i = 1:n,
        for j = 1:n,
            set(hh(i,j),'Visible','off');
        end;
    end;
end;


	