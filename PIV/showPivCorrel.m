function showPivCorrel(I,x,y,u,v,vecdata)

fig1 = gcf;

fig2 = findobj('Tag','CorrelFigure');
if (isempty(fig2)),
  pos1 = get(fig1,'Position');

  fig2 = figure('Position',[pos1(1)+pos1(3) pos1(2) pos1(3)/2 pos1(4)/2],...
                'Tag','CorrelFigure');
end;

figure(fig1);
set(fig1,'Renderer','painters');
data.imHandle = imshow(I,'n');
data.vecHandle = addquiverc(x,y,u,v,'y');

set(data.vecHandle,'ButtonDownFcn',@spcButtonDown);
set(data.imHandle,'ButtonDownFcn',@spcButtonDown);
set(fig1,'KeyPressFcn',@spcKeyPress);

data.fig1 = fig1;
data.ax1 = gca;
data.fig2 = fig2;

data.x = x;
data.y = y;
data.u = u;
data.v = v;

data.correl = vecdata.Correlation;
data.cx = (1:size(data.correl,1)) - ceil(size(data.correl,1)/2) - 1;
data.cy = (1:size(data.correl,2)) - ceil(size(data.correl,2)/2) - 1;

guidata(fig1,data);

uiwait(fig1);

set(fig1,'KeyPressFcn','');
set(data.vecHandle,'ButtonDownFcn','');

% ---------------------------------------------------
function spcButtonDown(obj,eventdata)

data = guidata(obj);

c = get(data.ax1,'CurrentPoint');

xd = [1;1]*data.x(:)' + [0;1]*data.u(:)';
yd = [1;1]*data.y(:)' + [0;1]*data.v(:)';

[m,ind] = min((xd(:)-c(1,1)).^2 + (yd(:)-c(1,2)).^2);
[q,ind] = ind2sub(size(xd),ind);
[i,j] = ind2sub(size(data.x),ind);

h = findobj(data.fig1,'Tag','vecHighlight');
if (~isempty(h)),
  delete(h);
end;
h = addplot(data.x(i,j),data.y(i,j),'ro','Tag','vecHighlight');

figure(data.fig2);
contourf(data.cx,data.cy,data.correl(:,:,i,j),20);
shading flat;
colorbar;
axis equal ij;

addquiverc(0,0,data.u(i,j),data.v(i,j),'k','AbsScale',1);

figure(data.fig1);

% ---------------------------------------------------
function spcKeyPress(obj,eventdata)

data = guidata(obj);

c = get(data.fig1,'CurrentCharacter');

switch lower(c),
 case {'q',char(13)},
  uiresume(data.fig1);
end;



