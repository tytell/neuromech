function varargout = spike2browser(varargin)

opt.margin = [0.12 0.02 0.1 0.02];
opt.gap = 0.01;
opt.showamptub = false;
opt.amptubchan = 'S1';
opt.showcursors = false;
opt.maxsamps = 5e5;

p = 1;
fig = -1;
if ((nargin >= 1) && (numel(varargin{1}) == 1) && ishandle(varargin{1}) && ...
        strcmp(get(varargin{1},'Type'),'figure'))
    fig = varargin{1};  
    p = 2;
end;

pn = [];
file = [];
F = {};
varnames = {};
if (nargin < p)
    file = [];
elseif (~ischar(varargin{p}))
    file = [];
    if (iscell(varargin{p}))
        F = varargin{p};
        varnames = cell(size(F));
        p = p+1;
    end;
elseif (strcmp(varargin{p},'base'))
    file = 'base';
    F = getvar('-all','-tostruct');
    varnames = fieldnames(F);
    F = struct2cell(F);
    p = p+1;
elseif (exist(varargin{p},'dir'))
    pn = varargin{p};
    file = [];
    p = p+1;
elseif (exist(varargin{p},'file'))
    file = varargin{p};
    p = p+1;
end;

opt = parsevarargin(opt,varargin(p:end),p+1,'typecheck',false);

if ((nargin == 0) || (isempty(file) && isempty(F)))
    filtspec = fullfile(pn,'*.mat');
    [fn,pn] = uigetfile(filtspec, 'Open Spike2 converted file');
    if (isempty(fn))
        fprintf('Cancelled...\n');
        return;
    end;
    file = fullfile(pn,fn);
end;

if (isempty(F))
    F = load(file);
    varnames = fieldnames(F);
    F = struct2cell(F);
end;

spike2chans = cellfun(@(x) (isstruct(x) && isfield(x,'title')), F);
F = F(spike2chans);
varnames = varnames(spike2chans);

names = cellfun(@(x) (x.title), F, 'UniformOutput',false);
good = ~ismember(names,{'Keyboard'});

F = F(good);
varnames = varnames(good);

if (islogical(opt.showamptub) && opt.showamptub)
    opt.showamptub = 'S1';
elseif (~ischar(opt.showamptub))
    opt.showamptub = '';
end;

if (~isempty(opt.showamptub))
    s1 = find(strcmp(names,opt.showamptub));
    
    if (~isempty(s1))
        [amp,tub,t] = extractAmpTub(F{s1});
        
        amp1 = F{s1};
        amp1.values = amp;
        amp1.title = 'ampullary';
       
        tub1 = F{s1};
        tub1.values = tub;
        tub1.title = 'tuberous';
        
        F = [F(1:s1); amp1; tub1; F(s1+1:end)];
        varnames = [varnames(1:s1); cell(2,1); varnames(s1+1:end)];
    end;
end;

n = length(F);

width = 1 - opt.margin(1) - opt.margin(2);
height = (1 - opt.margin(3)-opt.margin(4) - (n-1)*opt.gap)/n;
h = -1*ones(n,1);
hln = -1*ones(n,2);

if (opt.showcursors)
    tmark = F{1}.start + [0.33*F{1}.length 0.66*F{1}.length]*F{1}.interval;
end;

if (fig == -1)
    fig = figure;
else
    figure(fig);
end;
clf;
for i = 1:n,
    h(i) = axes('Position',[opt.margin(1) 1-opt.margin(4)-i*height-(i-1)*opt.gap width height], ...
        'Parent',fig);

    if (F{i}.length > opt.maxsamps)
        fac = floor(F{i}.length/opt.maxsamps);
    else
        fac = 1;
    end;
    
    if (isfield(F{i},'interval') && isfield(F{i},'values'))
        t1 = (0:fac:F{i}.length-1)*F{i}.interval;
        if (isfield(F{i},'start'))
            t1 = t1 + F{i}.start;
        end;
        plot(h(i), t1,F{i}.values(1:fac:end), 'HitTest','off');
        axis tight;
    elseif (isfield(F{i},'times'))
        tt = repmat(F{i}.times(:)',[3 1]);
        yy = repmat([0; 1; NaN],[1 F{i}.length]);
        plot(h(i), tt(:),yy(:), 'HitTest','off');
        ytick off;
        axis tight;
    end;
    
    varname1 = varnames{i};
    ind = last(varname1 == '_');
    if (~isempty(ind))
        varname1 = varname1(ind+1:end);
    end;
    if (isfield(F{i},'units'))
        ylabel(h(i),{varname1,F{i}.title, F{i}.units});
    else
        ylabel(h(i),{varname1,F{i}.title});
    end;
    
    if (i ~= n)
        xtick labeloff;
    else
        xlabel('Time (s)');
    end;
    
    if (opt.showcursors)
        hln(i,:) = vertplot(tmark,'k-');
    end;
end;

linkaxes(h,'x');
zoom xon;

if (~isempty(file))
    [~,fn] = fileparts(file);
    set(fig,'Name',fn);
end;

if (opt.showcursors)
    if (nargout == 0)
        set(fig, 'CloseRequestFcn',{@savebeforeclose,hln,file});
        
        uicontrol(fig, 'Style','pushbutton','String','Export', ...
            'Position',[20 20 100 35], 'Callback',{@savemarkers,hln,file});
    else
        set(fig, 'CloseRequestFcn',{@savebeforeclose,hln,''}, 'ToolBar','figure');
        
        oldunits = get(fig,'Units');
        set(fig,'Units','pixels');
        pos = get(fig,'Position');
        set(fig,'Units',oldunits);
        
        uicontrol(fig, 'Style','pushbutton','String','OK', ...
            'Position',[20 20 100 35], 'Callback',{@savemarkers,hln,''});
        uicontrol(fig, 'Style','pushbutton','String','Cancel', ...
            'Position',[pos(3)-20-100 20 100 35], 'Callback',@cancelbutton);
    end;
    
    set(hln,'ButtonDownFcn',{@markbuttondown,hln},'UserData',false);
    set(h, 'ButtonDownFcn',{@axbuttondown,hln});
    
    if (nargout > 0)
        set(fig,'UserData',false);
        uiwait(fig);
        
        if (all(ishandle(hln(:))) && get(fig,'UserData'))
            t0 = get(hln(1,1),'XData');
            t1 = get(hln(1,2),'XData');
            
            tmark = sort([t0(1) t1(1)]);
            varargout = {tmark};
        else
            varargout = {[]};
        end;
        delete(fig);
    end;
end;

function markbuttondown(src,evnt,hln)
% src - the object that is the source of the event
% evnt - empty for this property

set(src,'Selected','on');

[~,ind] = find(hln == src);

set(gcbf, 'WindowButtonMotionFcn', {@markdrag,src,ind,hln});
set(gcbf, 'WindowButtonUpFcn', @markbuttonup);



function markdrag(src,evnt, hln1,ind,hlnall)

ax = get(hln1, 'Parent');

pt = get(ax,'CurrentPoint');
set(hlnall(:,ind),'XData',[pt(1,1); pt(1,1)], 'UserData',true);



function markbuttonup(src,evnt)

set(gcbf, 'WindowButtonMotionFcn', []);
set(gcbf, 'WindowButtonUpFcn', []);


function axbuttondown(src,evnt,hln)
% src - the object that is the source of the event
% evnt - empty for this property

t0 = get(hln(1,1),'XData');
t1 = get(hln(1,2),'XData');

pt = get(src,'CurrentPoint');
if (abs(pt(1,1) - t0(1)) < abs(pt(1,1) - t1(1)))
    set(hln(:,1),'XData',[pt(1,1); pt(1,1)]);
else
    set(hln(:,2),'XData',[pt(1,1); pt(1,1)]);
end;


function savemarkers(src,evnt, hln, file)

t0 = get(hln(1,1),'XData');
t1 = get(hln(1,2),'XData');

tmark = sort([t0(1) t1(1)]);

if (~isempty(file))
    save(file,'tmark','-append');
else
    uiresume;
end;

set(hln,'UserData',false);
set(gcbf,'UserData',true);

function cancelbutton(src,evnt,dln)

set(gcbf,'UserData',false);
uiresume

function savebeforeclose(src,evnt, hln,file)
% User-defined close request function 
% to display a question dialog box 

if (~isempty(findobj(hln,'UserData',true)))
    selection = questdlg('Save markers before closing?',...
        'Save markers?',...
        'Yes','No','Cancel','Yes');
    switch selection,
        case 'Yes',
            savemarkers(src,evnt, hln,file);
            delete(gcf);
        case 'No',
            delete(gcf);
        case 'Cancel',
            % do nothing
    end
else
    delete(gcf);
end;



