function h = rose2(ang, varargin)

nbins = 24;
bins = linspace(0,2*pi,nbins+1);
opt.method = 'fill';
opt.opacity = 1;
opt.center = 0;
opt.join = false;
opt.normalize = false;
opt.tight = false;
opt.color = [];
opt.smooth = false;

if ((length(varargin) >= 1) && isnumeric(varargin{1})),
    if (numel(varargin{1}) == 1),
        nbins = varargin{1};
        bins = linspace(0,2*pi,nbins+1);
    else
        bins = varargin{1};
        if (bins(end) < 2*pi),
            bins(end+1) = 2*pi;
        end;
        nbins = length(bins);
    end;
    p = 2;
else
    p = 1;
end;

while (p <= length(varargin)),
    switch lower(varargin{p}),
      case {'line','fill'},
        opt.method = lower(varargin{p});
        p = p+1;
        
      case {'opacity','center','color'},
        opt.(lower(varargin{p})) = varargin{p+1};
        p = p+2;
        
      case {'join','normalize','tight','smooth'},
        opt.(lower(varargin{p})) = true;
        p = p+1;
        
      otherwise,
        error('Unrecognized option %s');
    end;
end;

ngroups = size(ang,2);
if (ngroups > 10),
    warning('rose2:largengroups',...
            'Number of groups is large. Did you forget to transpose?');
end;

n = zeros(length(bins),ngroups);

for i = 1:ngroups,
    n(:,i) = histc(mod(ang(:,i),2*pi),bins);
    
    if (opt.normalize),
        n(:,i) = n(:,i) / sum(n(:,i));
    end;
end;

n = n(1:end-1,:);

if (opt.join),
    if (opt.smooth),
        r = n;
        theta = (bins(1:end-1)+bins(2:end))/2;
        theta = theta(:);
        theta = theta(:,ones(1,ngroups));
    else,
        r = shiftdim(n,-1);
        r = r(ones(2,1),:,:);
        r = reshape(r,[2*nbins ngroups]);
        
        b1 = bins(1:end-1);
        b2 = bins([2:end-1 1]);
        theta = [b1(:)'; b2(:)'];
        theta = theta(:);
        theta = theta(:,ones(1,ngroups));
    end;
    
    if (isempty(opt.color)),
        col = 1:ngroups;
    else
        col = opt.color(:)';
    end;
else,
    r = zeros(4,nbins,ngroups);
    r([1 4],:,:) = opt.center;
    r(2,:,:) = shiftdim(n,-1);
    r(3,:,:) = shiftdim(n,-1);

    b1 = bins(1:end-1);
    b2 = bins([2:end-1 1]);
    theta = [b1(:)'; b1(:)'; b2(:)'; b2(:)'];
    theta = theta(:,:,ones(1,ngroups));
    
    if (isempty(opt.color)),
        col = shiftdim(1:ngroups,-1);
    else
        col = shiftdim(opt.color(:),-2);
    end;
    col = col(1,ones(1,nbins),:);
    
    r = reshape(r,[4 nbins*ngroups]);
    theta = reshape(theta,[4 nbins*ngroups]);
    col = reshape(col,[1 nbins*ngroups]);
end;

%deal with overwriting or not the plot
switch get(gca,'NextPlot'),
 case 'replace',
  cla reset;
  unhold = true;
  drawaxes = true;
 case 'replacechildren',
  cla;
  unhold = true;
  drawaxes = true;
 otherwise,
  unhold = false;
  drawaxes = false;
end;

%make a square plot briefly, just to see how Matlab autoscales in the
%range from 0 to lim in the space that will be available
if (drawaxes),
    lim = max(n(:)) + opt.center;
    curax = gca;
    axis equal off;
    pos = get(curax,'Position');
    testpos = [pos(1:2) min(pos(3:4))/2*[1 1]];
    axquick = axes('Position',testpos);
    plot([opt.center lim lim opt.center],[lim lim opt.center opt.center]);
    if (opt.tight),
        axis tight;
    end;
    set(axquick,'XTickLabel',{}, 'YTickLabel',{});
    xtick = get(axquick,'XTick');               % get the tick values
    lim = get(axquick,'XLim');
    lim = lim(2);

    delete(axquick);
    axes(curax);                            % return to the current axes

    %but save the tick positions
    tick = xtick(xtick >= 0);
    tick = tick - opt.center;

    %use a patch as the axis background
    circ = linspace(0,2*pi,51)';
    xcirc = cos(circ);
    ycirc = sin(circ);
    patch(lim*xcirc,lim*ycirc,'w');

    %draw the inner radius
    hold on;
    if (opt.center > 0),
        line(opt.center*xcirc, opt.center*ycirc,'Color','k');
    end;

    %draw the radial lines
    if (opt.center > 0),
        radial = 0:pi/6:2*pi;
        radial = radial(1:end-1);

        rr = [opt.center; lim];
        xradial = rr*cos(radial);
        yradial = rr*sin(radial);
    else
        radial = 0:pi/6:pi;
        radial = radial(1:end-1);

        rr = [-lim; lim];
        xradial = rr*cos(radial);
        yradial = rr*sin(radial);
    end;
    xradial(end+1,:) = NaN;
    yradial(end+1,:) = NaN;
    line(xradial(:),yradial(:),'Color','k','LineStyle',':');

    %and label them
    text(xradial(end-1,:)*1.1,yradial(end-1,:)*1.1,num2str(radial'*180/pi),...
         'HorizontalAlignment','center');
    if (opt.center == 0),
        text(xradial(1,:)*1.1,yradial(1,:)*1.1,num2str(radial'*180/pi+180),...
             'HorizontalAlignment','center');
    end;

    %and circular grid lines
    xcirc = repmat(xcirc,[1 length(tick)-2]) .* repmat(tick(2:end-1),[51 1]);
    ycirc = repmat(ycirc,[1 length(tick)-2]) .* repmat(tick(2:end-1),[51 1]);
    xcirc(end+1,:) = NaN;
    ycirc(end+1,:) = NaN;
    line(xcirc(:),ycirc(:),'Color','k','LineStyle',':');

    [q,circlabpos] = min(sum(n,2));
    circlabpos = sum(bins(circlabpos+[0 1]))/2;
    circlab = tick(2:end);
    circlab = num2str(circlab(:));
    text(tick(2:end)*cos(circlabpos), tick(2:end)*sin(circlabpos), circlab, ...
         'HorizontalAlignment','center','VerticalAlignment','middle');
end;

h = fill(r.*cos(theta), r.*sin(theta), col(ones(size(r,1),1),:));
set(h,'AlphaDataMapping','none', 'FaceAlpha',opt.opacity);

if (unhold),
    hold off;
end;

    
        