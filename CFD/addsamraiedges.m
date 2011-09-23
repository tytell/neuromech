function addsamraiedges(varargin)

opt.patchcolorder = 'kbrcm';
opt.keepaxislimits = true;
opt.showoverlap = false;

smallval = 1e-6;

if ((numel(varargin{1}) == 1) && ishandle(varargin{1}) && ...
        strcmp(get(varargin{1},'Type'),'axes')),
    [ax,V] = deal(varargin{1:2});
    p = 3;
else
    ax = gca;
    V = varargin{1};
    p = 2;
end;

opt = parsevarargin(opt,varargin(p:end),p, 'typecheck',false);

if (ischar(opt.patchcolorder) && (size(opt.patchcolorder,1) == 1)),
    opt.patchcolorder = opt.patchcolorder(:);
end;

if (opt.keepaxislimits),
    set(ax, 'XLimMode','manual','YLimMode','manual');
end;

xlo = cat(2,V.xlo);
ylo = xlo(2,:)';
xlo = xlo(1,:)';

xup = cat(2,V.xup);
yup = xup(2,:)';
xup = xup(1,:)';

patchlev = cat(1,V.level_number);

for i = 1:length(V),
    col = opt.patchcolorder(min(V(i).level_number+1,size(opt.patchcolorder,1)),:);
    
    xedge = cell(4,1);
    xedge{1} = [xlo(i); xup(i)];
    xedge{2} = [xup(i); xup(i)];
    xedge{3} = [xup(i); xlo(i)];
    xedge{4} = [xlo(i); xlo(i)];

    yedge = cell(4,1);
    yedge{1} = [ylo(i); ylo(i)];
    yedge{2} = [ylo(i); yup(i)];
    yedge{3} = [yup(i); yup(i)];
    yedge{4} = [yup(i); ylo(i)];
    
    if (opt.showoverlap),
        for j = [1 3],
            xe1 = xedge{j};
            ye1 = yedge{j};
            overlap = (patchlev == patchlev(i)+1) & ...
                (max(xe1) > xlo) & (min(xe1) < xup) & ...
                (ye1(1) >= ylo) & (ye1(1) <= yup);
            
            if (any(overlap)),
                if (xe1(end) > xe1(1))
                    sortmode = 'ascend';
                else
                    sortmode = 'descend';
                end;
                [xe1,ord] = sort([xlo(i); max(xlo(overlap),xlo(i)); ...
                    max(xlo(overlap),xlo(i))+smallval; min(xup(overlap),xup(i))-smallval; ...
                    min(xup(overlap),xup(i)); xup(i)], ...
                    sortmode);
                nover = sum(overlap);
                
                brks = [false; false(nover,1); true(nover,1); true(nover,1); ...
                    false(nover,1); false];
                brks = brks(ord);
                
                xe1(brks) = NaN;
                
                ye1 = ye1(1) * ones(size(xe1));
                
                xedge{j} = xe1;
                yedge{j} = ye1;
            end;
        end;
        
        for j = [2 4],
            xe1 = xedge{j};
            ye1 = yedge{j};
            overlap = (patchlev == patchlev(i)+1) & ...
                (max(ye1) > ylo) & (min(ye1) < yup) & ...
                (xe1(1) >= xlo) & (xe1(1) <= xup);
            
            if (any(overlap)),
                if (ye1(end) > ye1(1))
                    sortmode = 'ascend';
                else
                    sortmode = 'descend';
                end;
                [ye1,ord] = sort([ylo(i); max(ylo(overlap),ylo(i)); ...
                    max(ylo(overlap),ylo(i))+smallval; min(yup(overlap),yup(i))-smallval; ...
                    min(yup(overlap),yup(i)); yup(i)], ...
                    sortmode);
                nover = sum(overlap);
                
                brks = [false; false(nover,1); true(nover,1); true(nover,1); ...
                    false(nover,1); false];
                brks = brks(ord);
                
                ye1(brks) = NaN;
                
                xe1 = xe1(1) * ones(size(ye1));
                
                xedge{j} = xe1;
                yedge{j} = ye1;
            end;
        end;
    end;
    
    xedge1 = cat(1,xedge{:});
    yedge1 = cat(1,yedge{:});

    d = [1; abs(diff(xedge1)) + abs(diff(yedge1))];
    good = isnan(d) | (d > 2*smallval);
    
    xedge1 = xedge1(good);
    yedge1 = yedge1(good);
    
    good = [true; ~((isnan(xedge1(1:end-2)) | isnan(yedge1(1:end-2))) & ...
        (isnan(xedge1(3:end)) | isnan(yedge1(3:end)))); true];
    
    addplot(ax, xedge1(good),yedge1(good), 'Color',col);
end;
