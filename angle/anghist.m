function varargout = anghist(varargin)
% ANGHIST - Angular histogram plot
%
%          anghist(ang)
% or       anghist(ang,bins)
% or       anghist(angx,angy,bins)
%
% where angx and angy are cos and sin of the angle, respectively.  bins can
% be a scalar, in which case it determines the number of equally spaced bins
% (default 24).  If bins is a vector, it determines the edges of the bins,
% like in HIST
%
% SEE ALSO
%   HIST

% Mercurial revision hash: $Revision$ $Date$
% See http://rcn.ccs.tulane.edu/index.php5/Tytell_Matlab
% Copyright (c) 2012, Eric Tytell <tytell at jhu dot edu>

angx = [];
angy = [];

%look for the first option, if it exists
chararg = find(cellfun('isclass',varargin,'char'));
if (~isempty(chararg)),
    nnumeric = chararg(1)-1;
else
    nnumeric = length(varargin);
end;

if (nnumeric == 1),
    ang = varargin{1};
    bins = [];
    p = 2;
elseif (nnumeric == 2),
    if ((ndims(varargin{1}) == ndims(varargin{2})) && ...
            all(size(varargin{1}) == size(varargin{2}))),
        angx = varargin{1};
        angy = varargin{2};
        ang = atan2(angy,angx);
        bins = [];
    else
        ang = varargin{1};
        bins = varargin{2};
    end;
    p = 3;
elseif (nnumeric == 3),
    angx = varargin{1};
    angy = varargin{2};
    bins = varargin{3};

    ang = atan2(angy,angx);
    angx = nanmean(angx);
    angy = nanmean(angy);
    p = 4;
end;

%default values
%ctrsize = 0.75;
%bingap = 1/3;
ctrsize = 0.01;
bingap = 0.01;

%process options
while (p <= nargin),
    switch lower(varargin{p}),
     case 'ctrsize',
      ctrsize = varargin{p+1};
      p = p+2;
     case 'bingap',
      bingap = varargin{p+1};
      p = p+2;
     case 'freqlim',
      freqlim = varargin{p+1};
      p = p+2;

     otherwise,
      error('Unknown option %s.',varargin{p});
    end;
end;

if ((bingap < 0) || (bingap > 1)),
    error('Bin gap must be between 0 and 1.');
end;
if ((ctrsize <= 0) || (ctrsize >= 1)),
    error('Center size must be greater than 0 and less than 1.');
end;

%deal with overwriting or not the plot
switch get(gca,'NextPlot'),
 case 'replace',
  cla reset;
  unhold = 1;
 case 'replacechildren',
  cla;
  unhold = 1;
 otherwise,
  unhold = 0;
end;

ang = shiftdim(ang);
ang = mod(ang,2*pi);

if (isempty(bins)),
  bins = 24;
end;

%sort out the bin edges if we need to
if (numel(bins) == 1),
  edges = linspace(0,2*pi,bins+1)';
  db = (edges(2)-edges(1))/2;
  edges = edges-db;

  k = find(ang > 2*pi-db);
  ang(k) = ang(k)-2*pi;
else
  edges = shiftdim(bins);
end;

%define the bin centers 
ctrs = (edges(1:end-1) + edges(2:end))/2;
ctrs = repmat(ctrs,[1 size(ang,2)]);

%space out the different histogram elements from the bin center, so that
%they cover 1-bingap of the bin
db = mean(diff(edges));
spc = linspace(0,(1-bingap)*db,size(ang,2));
spc = spc - mean(spc);

ctrs = ctrs + repmat(spc,[size(ctrs,1) 1]);

ang(ang == 2*pi) = 0;

%actually do the histogram
n = histc(ang, edges);

%we shouldn't have anything in the 2*pi to Inf bin - if we do than
%something is wrong
if (any(n(end,:) > 0)),
  warning('Something weird is happening.');
end;
%get the total number and divide the number in each bin by it so we have
%a faction
n = n(1:end-1,:);
n = n./repmat(sum(n),[size(n,1) 1]);

% now set up the plot

%each "bar" is actually a thick line
%we define them in polar coordinates right now - a is the angle and r is
%the radius
a = repmat(shiftdim(ctrs,-1),[3 1 1]);
r = zeros(3,size(n,1),size(n,2));
r(2,:,:) = n;
r(3,:,:) = repmat(NaN,[1 size(n)]);

r = flatten(r,[1 2]);
a = flatten(a,[1 2]);

%ctrsize defines the fraction of the plot that is the center, but the
%true size depends on the maximum radius we'll be plotting.  c0 is the
%actual size of the center
lim = nanmax(r(:));
c0 = lim*ctrsize/(1-ctrsize);

%make a square plot briefly, just to see how Matlab autoscales in the
%range from 0 to lim in the space that will be available
curax = gca;
axis equal off;
pos = get(curax,'Position');
testpos = [pos(1:2) min(pos(3:4))*(1-ctrsize)/2*[1 1]];
axquick = axes('Position',testpos);
h = plot([0 lim lim 0],[lim lim 0 0]);
if (exist('freqlim','var')),
    axis([0 freqlim 0 freqlim]);
end;
xtick = get(axquick,'XTick');               % get the tick values
if (exist('freqlim','var') && (xtick(end) ~= freqlim)),
    xtick(end+1) = freqlim;
end;
delete(h);                              % delete the square plot
delete(axquick);
axes(curax);                            % return to the current axes

%but save the tick positions
tick = xtick(xtick >= 0);
tick = tick+c0;

%use a patch as the axis background
circ = linspace(0,2*pi,51)';
xcirc = cos(circ);
ycirc = sin(circ);
patch(tick(end)*xcirc,tick(end)*ycirc,'w');
hold on;
line(c0*xcirc, c0*ycirc,'Color','k');

radial = 0:pi/6:2*pi;
radial = radial(1:end-1);

rr = [c0 tick(end)]';
xradial = rr*cos(radial);
yradial = rr*sin(radial);
xradial(end+1,:) = NaN;
yradial(end+1,:) = NaN;
line(xradial(:),yradial(:),'Color','k','LineStyle','--');
text(xradial(end-1,:)*1.1,yradial(end-1,:)*1.1,num2str(radial'*180/pi),...
     'HorizontalAlignment','center');

xcirc = repmat(xcirc,[1 length(tick)-2]) .* repmat(tick(2:end-1),[51 1]);
ycirc = repmat(ycirc,[1 length(tick)-2]) .* repmat(tick(2:end-1),[51 1]);
xcirc(end+1,:) = NaN;
ycirc(end+1,:) = NaN;
line(xcirc(:),ycirc(:),'Color','k','LineStyle','--');

rnum = histc(ang(:),radial);
[q,rind] = min(rnum(1:end-1));
rlabel = radial(rind) + pi/12;
text(tick(2:end)*cos(rlabel),tick(2:end)*sin(rlabel), ...
     num2str(tick(2:end)'-c0));

if (isempty(angx)),
  angx = nanmean(cos(ang));
  angy = nanmean(sin(ang));
end;

h = plot((r+c0).*cos(a),(r+c0).*sin(a),'LineWidth',4);
numzero = ctrs;
for i = 1:size(n,2),
    good = n(:,i) > 0;
    numzero(good,i) = NaN;
end;
h2 = plot(c0.*cos(numzero),c0.*sin(numzero),'.');
h = cat(1,h,h2);

plot([zeros(size(angx)); angx]*c0,...
     [zeros(size(angy)); angy]*c0,'-');
plot([repmat(NaN,size(angx)); angx]*c0,...
     [repmat(NaN,size(angy)); angy]*c0,'o');

if (nargout == 1),
    varargout{1} = h;
end;

axis equal off;
if (unhold),
  hold off;
end;








