function summarizekinematics(gp,vel,varargin)

opt.varnames = [];

sz1 = size(varargin{1});

i = 1;
opt1 = 1;
while (i <= length(varargin)),
    if (isnumeric(varargin{i}) && (ndims(varargin{i}) == length(sz1)) && ...
        all(size(varargin{i}) == sz1)),
        opt1 = i+1;
        i = i+1;
    elseif (ischar(varargin{i})),
        switch lower(varargin{i}),
          case {'varnames'},
            opt.(lower(varargin{i})) = varargin{i+1};
            i = i+2;
            
          otherwise,
            error('Unrecognized option %s',varargin{i});
        end;
    end;
end;
        
vars = cat(3,varargin{1:opt1-1});
nvars = size(vars,3);

if (isempty(opt.varnames)),
    opt.varnames = cell(1,nvars);
    for i = 1:opt1-1,
        opt.varnames{i} = inputname(i+2);
    end;
end;

figure(1);

a = floor(sqrt(nvars));
b = ceil(sqrt(nvars));
for i = 1:nvars,
    subplot(a,b,i);
    plotgroups(vel,vars(:,:,i),gp,{'cmf'},'legend');
    title(opt.varnames{i});
end;

setnum = 1:size(vars,2);
setnum = setnum(ones(1,size(vars,1)),:);

X = cat(2,vel(:),reshape(vars,[size(vars,1)*size(vars,2) nvars]));
good = all(isfinite(X),2);

X = X(good,:);

[pc,scores,vnc,t2] = princomp(zscore(X));
flipsign = max(pc) < -min(pc);
pc(:,flipsign) = -pc(:,flipsign);
scores(:,flipsign) = -scores(:,flipsign);

figure(2);
subplot(3,3,1);
plot(vnc,'o-');
xlabel('PC');
ylabel('Variance');

subplot(3,3,4);
biplot(pc(:,1:2),'VarLabels',{'vel',opt.varnames{:}});

subplot(3,3,7);
biplot(pc(:,2:3),'VarLabels',{'vel',opt.varnames{:}});
xlabel('PC2');
ylabel('PC3');

subplot(2,3,2:3);
plotgroups(scores(:,1),scores(:,2),gp(good),{'cmf'},'legend');
xlabel('PC1');
ylabel('PC2');

subplot(2,3,5:6);
plotgroups(scores(:,2),scores(:,3),gp(good),{'cmf'},'legend');
xlabel('PC2');
ylabel('PC3');

figure(3);
normality = NaN(size(vel));
%use as the "normality" index a sum of the higher order PCs
%take everything with a variance component > 0.25
normpc = vnc > 0.25;
%never include PC1 and always include PC2
normpc([1 2]) = [false true];

normality(good) = sqrt(sum(scores(:,normpc).^2, 2));

a = floor(sqrt(nvars));
b = ceil(sqrt(nvars));
for i = 1:nvars,
    subplot(a,b,i);
    scatter(vel(:),flatten(vars(:,:,i)),pi*4^2,normality(:),'fill');
    title(opt.varnames{i});
end;


