function pivcorrplot(x,y,u,v,corr)

data.x = x;
data.y = y;
data.u = u;
data.v = v;
data.u1 = u(:,:,1);
data.v1 = v(:,:,1);
data.col1 = 'k';
data.corr = corr;
data.corr1 = data.corr(:,:,1,1,1);
data.sel = [1 1];

data.DU = velderiv(x,y,u,v);
data.bg = cat(3,data.DU.dvdx) - cat(3,data.DU.dudy);
data.bg1 = data.bg(:,:,1);
data.bgx = x(1,:);
data.bgy = y(:,1);

data.Figure = findobj('Tag','pivplotMain');
if (isempty(data.Figure)),
    data.Figure = figure('DoubleBuffer','on','Renderer','painters',...
                         'Tag','pivplotMain',...
                         'KeyPressFcn',@KeyPress);
    data.Axes = axes('Parent',data.Figure);
else
    data.Axes = findobj(data.Figure,'Type','axes');
end;

data.corrFigure = findobj('Tag','pivplotCorr');
if (isempty(data.corrFigure)),
    pos = get(data.Figure,'Position');
    data.corrFigure = figure('Position',[pos(1)+pos(3) pos(2)+0.5*pos(4) ...
                        0.5*pos(3) 0.5*pos(4)],...
                             'DoubleBuffer','on','Tag','pivplotCorr',...
                             'KeyPressFcn',@KeyPress);
    data.corrAxes = axes('Parent',data.corrFigure);
else
    data.corrAxes = findobj(data.Figure,'Type','axes');
end;

data.hBg = -1;
data.hVec = -1;
data.hCorr = -1;

data = Update(data, [1 2 10]);
guidata(data.Figure, data);

uiwait(data.Figure);

set(data.Figure, 'KeyPressFcn',[]);
set(get(data.Axes,'Children'), 'ButtonDownFcn', [], ...
                  'UIContextMenu', []);

% ----------------------------------------
function data = Update(data, what)

if (any(what == 1)),                    % update background
    if (ishandle(data.hBg)),
        sz = [length(get(data.hBg,'YData')) ...
              length(get(data.hBg,'XData'))];
        if (all(sz == [size(data.bg,1) size(data.bg,2)])),
            set(data.hBg, 'CData', data.bg1);
            redraw = 0;
        else
            delete(data.hBg);
            redraw = 1;
        end;
    else
        redraw = 1;
    end;

    if (redraw),
        if (~isempty(data.bg1)),
            figure(data.Figure);

            hold on;
            if (isa(data.bg1,'uint8') | isa(data.bg1,'uint16')),
                data.hBg = imshow(data.bgx,data.bgy,data.bg1);
            else
                data.hBg = imagesc(data.bgx,data.bgy,data.bg1);
            end;
            hold off;

            ch = get(data.Axes,'Children');
            set(data.Axes,'Children',ch([end 1:end-1]));
        end;
    end;

    axis equal ij tight;
end;
if (any(what == 2)),                    % update vectors
    figure(data.Figure);
    if (ishandle(data.hVec)),
        delete(data.hVec);
    end;
    hold on;
    data.hVec = quiverc(data.x,data.y,data.u1,data.v1,data.col1);
    hold off;
    axis equal ij tight;

    set(data.hVec, 'ButtonDownFcn', @vecButtonDown);
end;
if (any(what == 3)),
    figure(data.Figure);
    if (ishandle(data.hSelVec)),
        delete(data.hSelVec);
    end;
    hold on;
    data.hSelVec = 
if (any(what == 10)),                   % update correlation window
    if (ishandle(data.hCorr)),
        set(data.hCorr,'CData',data.corr1);
    else
        figure(data.corrFigure);

        data.hCorr = imagesc(data.corr1);
        caxis([-1 1]);
        colormap default;
        colorbar;
    end;
end;

% ----------------------------------------
function KeyPress(hObj, eventdata)

data = guidata(hObj);

c = get(data.Figure, 'CurrentCharacter');

switch lower(c),
 case 'q',
  uiresume(data.Figure);
end;

% ----------------------------------------
function vecButtonDown(hObj, eventdata)

data = guidata(hObj);

c = get(data.Axes, 'CurrentPoint');
[q,i] = min((c(1,1)-data.x(:)).^2 + (c(1,2)-data.y(:)).^2);

[i,j] = ind2sub(size(data.x),i);

data.sel = [i j];
data.corr1 = data.corr(:,:,i,j,1);

data = Update(data, 10);
guidata(hObj, data);




