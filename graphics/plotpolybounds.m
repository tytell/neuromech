function h = plotpolybounds(x,y,varargin)
% function plotPolyBounds(x,y,[alpha])
% Plots a regression line and (1-alpha)% confidence interval on the current
% figure.  Alpha is optional; by default it's 0.05.

opt.showstats = true;

if ((nargin >= 3) && isnumeric(varargin{1}) && (numel(varargin{1}) == 1)),
    alpha = varargin{1};
    p = 2;
else
    alpha = 0.05;
    p = 1;
end;

opt = parsevarargin(opt,varargin(p:end),p);

% Number of point to calculate for the confidence interval
n = 20;

% make sure we have column vectors
x = shiftdim(x);
y = shiftdim(y);

% Get rid of NaNs
k = find(isfinite(x) & isfinite(y));
x = x(k);
y = y(k);

% Do the regression in two ways -- the first one gives us r2 and P values,
% and the second gives us the S variable that we use to calculate the
% bounds
[~,bint,~,~,stats] = regress(y,[x ones(size(x))]);
[b,S] = polyfit(x,y,1);
sem = (bint(:,2) - bint(:,1))/2 / tinv(1-alpha,S.df);

xx = linspace(min(x),max(x),n);
[yy,dy] = polybounds(xx,b,S,alpha,1,1);

% plot the lines
hold on;
h(1) = plot(xx,yy,'k-','LineWidth',2);
h(2) = plot([xx xx(end:-1:1)],[yy+dy yy(end:-1:1)-dy],'-','Color',[0.8 0.8 0.8]);
hold off;

if (opt.showstats),
    % Show P and r2 values
    text(0.2,0.8,sprintf('P = %.3f',stats(3)),'Units','normalized');
    text(0.2,0.7,sprintf('r2 = %.3f',stats(1)),'Units','normalized');
end;

text(0.2,0.6,sprintf('y = (%.3g+-%.3g)x + (%.3g+-%.3g)',b(1),sem(1),b(2),sem(2)),'Units','normalized');
