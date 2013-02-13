function podkinematics(t,x,y,width,varargin)

opt.smoothcurve = 0;
opt.debug = true;
opt.ncomponents = 7;
opt.groups = [];
opt.col = 'krgbcmy';

comx = [];
comy = [];

if ((size(width,1) == 1) && (size(width,2) == size(x,2)) && ...
    size(width,3) == size(x,3) && (nargin >= 5) && ...
    isnumeric(varargin{1}) && (ndims(width) == ndims(varargin{1})) && ...
    all(size(width) == size(varargin{1}))),
    comx = width;
    comy = varargin{1};
    p = 2;
end;

if ((nargin >= 4+p) && isnumeric(varargin{p})),
    len = varargin{p};
    p = p+1;
else
    len = [];
    p = 1;
end;

opt = parsevarargin(opt,varargin(p:end),'firstoptionnumber',4+p,'typecheck',false);

nstartvars = 3;
npt = size(x,1);
nfr = size(x,2);
nreps = size(x,3);

if (numel(len) == nreps),
    len = reshape(len,[1 1 nreps]);
    
    %normalize everything by body length, if it's passed in
    x = x ./ len(ones(npt,1),ones(nfr,1),:);
    y = y ./ len(ones(npt,1),ones(nfr,1),:);
    
    if (~isempty(comx)),
        comx = comx ./ len(1,ones(nfr,1),:);
        comy = comy ./ len(1,ones(nfr,1),:);
    end;
end;

if (isempty(comx)),
    [comx,comy] = ctrofmasspos(x,y,width);
end;

comvx = deriv(t,comx,2);
comvy = deriv(t,comy,2);

comspeed = sqrt(comvx.^2 + comvy.^2);

[s,bodycurve,sp] = curvature(x,y,1,'spline','smooth',opt.smoothcurve);
bodycurvemn = nanmean(bodycurve(:));
bodycurvestd = nanstd(bodycurve(:));

fr = 1:nfr*nreps;
good = all(isfinite(bodycurve));

dxy = NaN(2,nfr,nreps);
dxy(:,good) = fnval(fnder(sp,[1 0]),{0,fr(good)});

headang = atan2(dxy(2,:,:),dxy(1,:,:));
headang(isnan(comspeed)) = NaN;
headang = unwrap(headang);
headang = headang - repmat(nanmean(headang,2),[1 size(headang,2) 1]);

h = nanmean(sqrt(diff(comx).^2 + diff(comy).^2));

pathcurve = curvature(comx(1,:,:),comy(1,:,:),2,'splineindiv',...
                             'smooth',h(1,:,:)/5);

bodycurven = (bodycurve - bodycurvemn) / bodycurvestd;

X = cat(1,comspeed,pathcurve,headang,bodycurven);
sz = size(X);
X = reshape(X,[sz(1) sz(2)*sz(3)]);
X = X';

Xmean1 = nanmean(X(:,1:nstartvars));
Xstd1 = nanstd(X(:,1:nstartvars));

X(:,1:nstartvars) = (X(:,1:nstartvars) - Xmean1(ones(nfr*nreps,1),:)) ./ ...
    Xstd1(ones(nfr*nreps,1),:);

good = all(isfinite(X),2);

varnames = cell(1,size(X,2));
varnames{1} = 'speed';
varnames{2} = 'pathcurve';
varnames{3} = 'headang';
for i = 4:length(varnames),
    varnames{i} = sprintf('curve%d',i-2);
end;

[coefs,scores,vars,t2] = princomp(X(good,:));
scores = NaN(size(X));
scores(good,:) = X(good,:) * coefs;

scores = scores';
scores = reshape(scores,sz);

if (opt.debug),
    h = -1*ones(opt.ncomponents,3);
    lo = Inf(1,3);
    hi = -Inf(1,3);
    
    subplot(1,8,1);
    barh(vars(1:opt.ncomponents));
    axis ij tight;

    u = unique(opt.groups(~isnan(opt.groups)));
    
    for i = 1:opt.ncomponents,
        h(i,1) = subplot(opt.ncomponents,6,6*(i-1) + 2);
        bar(coefs(1:nstartvars,i)');
        axis tight;
        set(gca,'Box','off','XTickLabel',[]);
        lo(1) = min([lo(1); coefs(1:nstartvars,i)]);
        hi(1) = max([hi(1); coefs(1:nstartvars,i)]);
        
        h(i,2) = subplot(opt.ncomponents,6,6*(i-1) + (3:4));
        hln = plot(1:size(X,2)-nstartvars,coefs(nstartvars+1:end,i)','.-');
        axis tight;
        set(hln,'Clipping','off');
        set(gca,'XTick',[]);
        lo(2) = min([lo(2); coefs(nstartvars+1:end,i)]);
        hi(2) = max([hi(2); coefs(nstartvars+1:end,i)]);

        h(i,3) = subplot(opt.ncomponents,6,6*(i-1) + (5:6));
        good = any(isfinite(scores(i,:,:)),3);
        if (~isempty(opt.groups)),
            cla;
            hold on;
            
            for j = 1:length(u),
                isgp = opt.groups == u(j);
                plot(squeeze(t(1,good,isgp)),squeeze(scores(i,good,isgp)),opt.col(j));
            end;
        else
            plot(squeeze(t(1,good,:)),squeeze(scores(i,good,:)));
        end;
        
        axis tight;
        set(gca,'XTickLabel',[]);
        lo(3) = min([lo(3); flatten(scores(i,good,:))]);
        hi(3) = max([hi(3); flatten(scores(i,good,:))]);
    end;
    
    linkaxes(h(:,1),'xy');
    linkaxes(h(:,2),'xy');
    linkaxes(h(:,3),'xy');
    
    axes(h(end,1));
    ylim([lo(1) hi(1)]);
    set(gca,'XTickLabel',varnames(1:nstartvars));
    
    axes(h(end,2));
    ylim([lo(2) hi(2)]);
    xlabel('Arc length');
    
    axes(h(end,3));
    ylim([lo(3) hi(3)]);
    set(gca,'XTickMode','auto');
    xlabel('Time');
end;

        
        