function angboxplot(varargin)

color = {};
alpha = 0.05;
group = [];
if (nargin == 1),
    ang = varargin{1};
    x = repmat(1:size(ang,2),[size(ang,1) 1]);
elseif ((nargin == 2) & (numel(varargin{2}) == 1)),
    ang = varargin{1};
    x = repmat(1:size(ang,2),[size(ang,1) 1]);
    alpha = varargin{2};
else,
    x = varargin{1};
    ang = varargin{2};

    if (any(size(x) == 1) & (length(x) == size(ang,2))),
        x = repmat(shiftdim(x)',[size(ang,1) 1]);
    end;

    if (nargin >= 3),
        for i = 3:nargin,
            if (isnumeric(varargin{i}) & (numel(varargin{i}) == 1)),
                alpha = varargin{i};
            elseif (ischar(varargin{i}) | (numel(varargin{i}) == 3)),
                color = varargin{i};
            elseif (all(size(varargin{i}) == size(ang))),
                group = varargin{i};
            end;
        end;
    end;
end;

group = repmat(1:size(ang,2),[size(ang,1) 1]);

nonan = find(isfinite(x) & isfinite(ang));
x = x(nonan);
ang = ang(nonan);
group = group(nonan);

[gp,q,gpind] = unique([x(:) group(:)],'rows');
gpx = gp(:,1);

if (~ishold),
    cla reset;

    circ = linspace(0,2*pi,51)';
    patch(gpx(end)*cos(circ),gpx(end)*sin(circ),'w');
    line(gpx(1)*cos(circ),gpx(1)*sin(circ),'Color','k','LineStyle','-');

    hold on;
    line(cos(circ)*gpx(2:end-1)',sin(circ)*gpx(2:end-1)',...
         'Color','k','LineStyle','--');

    radial = 0:pi/6:2*pi;
    radial = radial(1:end-1);
    line(gpx([1 end])*cos(radial), gpx([1 end])*sin(radial),...
         'Color','k','LineStyle','--');
    text(gpx(end)*cos(radial)*1.1, gpx(end)*sin(radial)*1.1,...
         num2str(radial'*180/pi),'HorizontalAlignment','center');

    unhold = 1;
else
    unhold = 0;
end;

for i = 1:length(gpx),
    k = find(gpind == i);

    [a1,r1] = angmean(ang(k));
    [sem1,conf1] = angconf(ang(k),alpha);
    x1 = x(k(1));

    theta = conf1(1):pi/36:conf1(2);
    theta(end+1) = conf1(2);

    meanplot{1,i} = x1*cos(a1);
    meanplot{2,i} = x1*sin(a1);
    rx(:,i) = [0; r1*cos(a1)]*gpx(1);
    ry(:,i) = [0; r1*sin(a1)]*gpx(1);
    stdplot{1,i} = x1*cos(theta);
    stdplot{2,i} = x1*sin(theta);
    angplot{1,i} = x1*cos(ang(k));
    angplot{2,i} = x1*sin(ang(k));
end;

if (~isempty(color)),
    color = {'Color',color};
end;

% plot(angplot{:},'Marker','.','LineStyle','none');
plot(stdplot{:},'LineWidth',2,color{:});
plot(rx,ry,'-',color{:});
plot(meanplot{:},'Marker','o','LineStyle','none',color{:});

axis equal off;

if (unhold),
    hold off;
end;

