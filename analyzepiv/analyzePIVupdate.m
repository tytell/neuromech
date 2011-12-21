function data = analyzePIVupdate(data, what, options)

if (nargin == 2),
    options = {};
elseif (ischar(options)),
    options = {options};
end;

if (strcmp(what,'all')),
    what = {'vectors','regions','background'};
end;

if (ischar(what)),
    what = {what};
end;

charOpts = find(cellfun('isclass',options,'char'));
if (strmatch('rescale',options(charOpts))),
    data.vectorOpts.AbsScale = [];
end;

axes(data.Axes);
keephold = ishold;
hold on;

for i = 1:length(what),
    switch what{i},
        case 'vectors',
            delete(data.hVectors);
            if (size(data.PIV.x,3) > 1),
                data.vectorOpts.x = data.PIV.x(:,:,data.curFrame);
                data.vectorOpts.y = data.PIV.y(:,:,data.curFrame);
            else
                data.vectorOpts.x = data.PIV.x;
                data.vectorOpts.y = data.PIV.y;
            end;
            data.vectorOpts.u = data.PIV.u(:,:,data.curFrame) - ...
                data.subtractVector(1);
            data.vectorOpts.v = data.PIV.v(:,:,data.curFrame) - ...
                data.subtractVector(2);
            [data.hVectors,data.vectorOpts] = quiverc(data.vectorOpts);
            set(data.hVectors,'HitTest','off');

        case 'regions',
            %       if (~isempty(strnummatch('recalc',options))),
            %           rgn = strnummatch('rgn',options);
            %           if (~isempty(rgn)),
            %               rgn = options{rgn+1};
            %           end;
            %           data = feval(data.calcGuiFcns.apCalculate,data,[],[]);
            %       end;


        case 'background',
            if (isempty(data.Background)),
                if (ishandle(data.hBackground)),
                    delete(data.hBackground);
                end;
            else
                bg = data.PIV.(data.Background);
                if (size(bg,3) > 1),
                    bg = bg(:,:,data.curFrame);
                end;
                if (ishandle(data.hBackground)),
                    sz = size(get(data.hBackground,'CData'));
                    if (all(size(bg) == sz)),
                        set(data.hBackground,'CData',bg);
                    else
                        delete(data.hBackground);
                    end;
                end;
                % we use this weird if structure, rather than a normal else,
                % because we might potentially delete the background handle
                % in the above if
                if (~ishandle(data.hBackground)),
                    data.hBackground = imagesc(data.PIV.x(1,:,1),...
                        data.PIV.y(:,1,1),bg);
                    caxis auto;
                    if (data.bgSymmetric),
                        ca = caxis;
                        ca = max(abs(ca));
                        ca = [-ca ca];
                        caxis(ca);
                    end;
                    colorbar;
                end;
            end;
    end;
end;
if (~keephold),
    hold off;
end;

% set the stacking so that the background is the deepest, the vectors are
% above the background, and everything else is above the vectors
children = get(data.Axes,'Children');
indvec(1) = find(children == data.hVectors(1));
if (length(data.hVectors) == 2),
    indvec(2) = find(children == data.hVectors(2));
end;
indbg = find(children == data.hBackground);
indother = setdiff(1:length(children),[indbg indvec]);
set(data.Axes,'Children',children([indother indvec indbg]));

if (strmatch('rescale',options(charOpts))),
    set(data.ScaleEdit,'String',num2str(data.vectorOpts.Scale));
end;


% -------------------------------------------------
function data = apUpdateRgnVectors(data, rgns, recalc)

if ((nargin == 1) | isempty(rgns)),
    rgns = 1:size(data.Regions,1);
end;

if (size(data.PIV.x,3) > 1),
    px = data.PIV.x(:,:,data.curFrame);
    py = data.PIV.y(:,:,data.curFrame);
else
    px = data.PIV.x(:,:,1);
    py = data.PIV.y(:,:,1);
end;

d = nanmean([flatten(diff(px,[],2)); ...
             flatten(diff(py,[],1))]);

keephold = ishold;
hold on;
for i = 1:length(rgns),
    rgn = rgns(i);

    if (recalc),
        h = data.Regions(rgn,data.curFrame).handle;

        x0 = get(h,'XData');
        y0 = get(h,'YData');

        s0 = [0 cumsum(sqrt(diff(x0).^2 + diff(y0).^2))];

        pu = data.PIV.u(:,:,data.curFrame);
        pv = data.PIV.v(:,:,data.curFrame);
        if data.PIV.isStereo
            pw = data.PIV.w(:,:,data.curFrame);
        else
            pw = [];
        end

        s1 = 0:d:s0(end);
        x1 = interp1(s0,x0, s1);
        y1 = interp1(s0,y0, s1);

        u1 = interp2(px,py,pu, x1,y1);
        v1 = interp2(px,py,pv, x1,y1);
        w1 = interp2(px,py,pw, x1,y1); % [] if pw = []

        data.Regions(rgn,data.curFrame).xcontour = x1;
        data.Regions(rgn,data.curFrame).ycontour = y1;
        data.Regions(rgn,data.curFrame).ucontour = u1;
        data.Regions(rgn,data.curFrame).vcontour = v1;
        data.Regions(rgn,data.curFrame).wcontour = w1;
    end;

%     if (ishandle(data.Regions(rgn,data.curFrame).hVectors)),
%         delete(data.Regions(rgn,data.curFrame).hVectors);
%     end;
%     
%     opts = data.vectorOpts;
%     opts.AbsScale = opts.Scale;
%     opts.Show = sqrt(opts.Show);
% 
%     if (rgn == data.curRgn),
%         opts.Color = 'r';
%     else
%         opts.Color = 'k';
%     end;
%     opts.x = data.Regions(rgn,data.curFrame).xcontour;
%     opts.y = data.Regions(rgn,data.curFrame).ycontour;
%     opts.u = data.Regions(rgn,data.curFrame).ucontour;
%     opts.v = data.Regions(rgn,data.curFrame).vcontour;
%     data.Regions(rgn,data.curFrame).hVectors = quiverc(opts);
%     set(data.Regions(rgn,data.curFrame).hVectors, 'HitTest','off');
end;

if (recalc),
%    data = feval(data.calcGuiFcns.apCalculate,data,[],[]);
end;

if (~keephold),
    hold off;
end;
