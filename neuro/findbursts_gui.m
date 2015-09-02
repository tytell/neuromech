function data = findbursts_gui(data, varargin)

opt.interburstdur = [];
opt.threshold = [];
opt.minspikes = [];
opt.goodchan = [];
opt.override = struct([]);
opt.show = 5;
opt.quiet = false;
opt.eventtimes = [];

if ((nargin >= 2) && isnumeric(data) && isnumeric(varargin{1}) && ...
        any(size(varargin{1}) == length(data)))
    data = struct('t',data, 'sig',varargin{1});
    args = varargin(2:end);
    p = 2;
    isdata = false;
else
    args = varargin;
    p = 1;
    isdata = true;
end
opt = parsevarargin(opt, args, p+1, 'typecheck',false);

if (~opt.quiet)
    fig = openfig(mfilename, 'new');
    set(fig,'WindowStyle','normal');
    gdata = guihandles(fig);
    set(gdata.prevButton,'Enable','off');
else
    gdata = struct;
end

nchan = size(data.sig,2);
gdata.chan = 1;
if (size(opt.threshold,1) == 2) && (size(opt.threshold,2) == nchan)
    gdata.thresh = opt.threshold;
elseif length(opt.threshold) == nchan
    gdata.thresh = [-1; 1]*abs(opt.threshold(:)');
else
    gdata.thresh = [-1; 1] * 2*nanstd(data.sig);
end
if (isempty(opt.interburstdur) || (length(opt.interburstdur) ~= nchan))
    gdata.interburst = 0.3*ones(1,nchan);
else
    gdata.interburst = opt.interburstdur;
end
if (isempty(opt.minspikes) || (length(opt.minspikes) ~= nchan))
    gdata.minspikes = 2*ones(1,nchan);
else
    gdata.minspikes = opt.minspikes;
end

if (isempty(opt.goodchan) || (length(opt.goodchan) ~= nchan))
    data.goodchan = true(1,nchan);
else
    data.goodchan = opt.goodchan > 0;
end

if (isfield(data,'burst'))
    data = rmfield(data,'burst');
end
gdata.data = data;
gdata.hspikes = -1;
gdata.hbursts = -1;
gdata.hdiag = [-1; -1];
gdata.hburstrate = -1;
gdata.hburstdur = -1;
gdata.selburstind = [];
gdata.hselburst = [-1 -1 -1 -1 -1];
gdata.burstoverride = repmat(struct('on',[],'off',[],'remove',[]),[1 nchan]);
gdata.data.spiket = cell(1,nchan);
gdata.data.spikeamp = cell(1,nchan);
gdata.eventtimes = opt.eventtimes;

gdata = update_spikes(gdata, true);

if ~isempty(opt.override)
    good = true;
    for c = 1:nchan
        if length(opt.override(c).on) < length(gdata.data.burst(c).on)
            good = false;
            break;
        end
    end
    if ~good
        warning('Burst override parameter is different from the number of bursts detected. Ignoring...');
    else
        gdata.burstoverride = opt.override;
    end
end
if (~opt.quiet)
    setUICallbacks(gdata);
    
    gdata = show_plot(gdata, 'all');
    guidata(fig,gdata);
    
    uiwait(fig);
    
    gdata = guidata(fig);
    if (ishandle(fig))
        delete(fig);
    end
    f = findobj('Tag','burstrate');
    if (~isempty(f))
        delete(f);
    end
end

data = gdata.data;
data.burst = get_burst_override(gdata,1:nchan);
for i = 1:nchan
    [c,ord] = sort(data.burst(i).ctr);
    data.burst(i).on = data.burst(i).on(ord);
    data.burst(i).off = data.burst(i).off(ord);
    data.burst(i).ctr = data.burst(i).ctr(ord);
    data.burst(i).nspike = data.burst(i).nspike(ord);
    data.burst(i).isover = data.burst(i).isover(ord);
    
    good = ~isnan(data.burst(i).ctr);
    data.burst(i).on = data.burst(i).on(good);
    data.burst(i).off = data.burst(i).off(good);
    data.burst(i).ctr = data.burst(i).ctr(good);
    data.burst(i).nspike = data.burst(i).nspike(good);
    data.burst(i).isover = data.burst(i).isover(good);
end
if any(cat(2,data.burst.isover))
    data.override = gdata.burstoverride;
else
    data.override = struct([]);
end

data.spikethreshold = gdata.thresh;
data.interburstdur = gdata.interburst;
data.minspikes = gdata.minspikes;

for i = 1:length(data.spiket)
    if isempty(data.spiket{i})
        data.spiket{i} = NaN;
        data.spikeamp{i} = NaN;
    end
end

spiket1 = catuneven(2,data.spiket{data.goodchan});
data.spiket = NaN(size(spiket1,1),nchan);
data.spiket(:,data.goodchan) = spiket1;
spikeamp1 = catuneven(2,data.spikeamp{data.goodchan});
data.spikeamp = NaN(size(spiket1,1),nchan);
data.spikeamp(:,data.goodchan) = spikeamp1;

for i = 1:length(data.burst)
    if isempty(data.burst(i).ctr)
        data.burst(i).ctr = NaN;
        data.burst(i).on = NaN;
        data.burst(i).off = NaN;
    end
    if isempty(data.burst(i).nspike)
        data.burst(i).nspike = NaN;
    end
end
ctr1 = {data.burst(data.goodchan).ctr};
ctr1 = cellfun(@(x) x', ctr1, 'UniformOutput',false);
burstt1 = catuneven(2,ctr1{:});
data.burstt = NaN(size(burstt1,1),nchan);
data.burstt(:,data.goodchan) = burstt1;

on1 = {data.burst(data.goodchan).on};
on1 = cellfun(@(x) x', on1, 'UniformOutput',false);
data.burston = NaN(size(data.burstt));
data.burston(:,data.goodchan) = catuneven(2,on1{:});

off1 = {data.burst(data.goodchan).off};
off1 = cellfun(@(x) x', off1, 'UniformOutput',false);
data.burstoff = NaN(size(data.burstt));
data.burstoff(:,data.goodchan) = catuneven(2,off1{:});

if isfield(data,'phase') && (~isfield(data,'amp') || (data.amp > 0))
    data.spikephase = NaN(size(data.spiket));
    data.burstphase = NaN(size(data.burstt));
    data.spikecyclet = NaN(size(data.spiket));
    data.burstcyclet = NaN(size(data.burstt));
    data.burstcycle = NaN(size(data.burstt));
    data.burststimphase = NaN(size(data.burstt));
    data.burststimcycle = NaN(size(data.burstt));

    uphase = unwrap(2*pi*data.phase) / (2*pi);
    goodphase = isfinite(uphase) & [true; diff(uphase) > 0];
    
    for i = 1:size(data.spiket,2)
        if (data.goodchan(i))
            good = isfinite(data.spiket(:,i));
            if any(good)
                data.spikephase(good,i) = interp1(data.t(goodphase),uphase(goodphase), data.spiket(good,i));
                good = isfinite(data.spikephase(:,i));

                cycle1 = floor(data.spikephase(good,i));
                data.spikecyclet(good,i) = interp1(uphase(goodphase),data.t(goodphase), cycle1);
            end

            good = isfinite(data.burstt(:,i));
            if any(good)
                data.burstphase(good,i) = interp1(data.t(goodphase),uphase(goodphase), data.burstt(good,i));
                good = isfinite(data.burstphase(:,i));

                if isfield(data,'stimphase')
                    goodstimphase = isfinite(data.stimphase) & [true; diff(data.stimphase) > 0];
                    data.burststimphase(good,i) = interp1(data.t(goodstimphase),data.stimphase(goodstimphase), data.burstt(good,i));
                    data.burststimcycle(good,i) = interp1(data.t(goodstimphase),data.stimcycle(goodstimphase), data.burstt(good,i));
                end
            end
        end
    end

    data.spikephase = mod(data.spikephase,1);
    data.burstphase = mod(data.burstphase,1);
end

if isfield(data,'stimfreq')
    if (length(data.stimfreq) == 1)
        ncycle = max(data.cycle);
        data.stimcyclet = (0:1/data.stimfreq:ncycle-1)';
    else
        data.stimcyclet = interp1(uphase(goodphase),data.t(goodphase),floor(min(uphase)):ceil(max(uphase)), ...
            'linear','extrap')';
    end

    [data.burstspercycle,data.burstcycle] = histc(data.burstt,data.stimcyclet);
    data.burstcycle(data.burstcycle == 0) = NaN;
    good = isfinite(data.burstcycle);
    data.burstcyclet = NaN(size(data.burstcycle));
    data.burstcyclet(good) = data.stimcyclet(data.burstcycle(good));
end

if (~opt.quiet)
    ibdtxt = sprintf('%g ',gdata.interburst);
    thtxt = sprintf('%g ',gdata.thresh(1,:));
    thtxt = [thtxt(1:end-1) ';' sprintf('%g ',gdata.thresh(2,:))];
    mstxt = sprintf('%g ',gdata.minspikes);
    goodtxt = sprintf('%d ',gdata.data.goodchan);
    
    if isdata
        fprintf(['%s = findbursts_gui(%s, ''threshold'', [%s], ''interburstdur'', [%s],' ...
            '''minspikes'', [%s], ''goodchan'', [%s], ''quiet'')\n'], ...
            inputname(1), inputname(1), thtxt(1:end-1), ibdtxt(1:end-1), mstxt(1:end-1), goodtxt(1:end-1));
    else
        fprintf(['data = findbursts_gui(%s,%s, ''threshold'', [%s], ' ...
            '''interburstdur'', [%s], ''minspikes'', [%s], ''goodchan'', [%s], ''quiet'')\n'], ...
            inputname(1), inputname(2),...
            thtxt(1:end-1), ibdtxt(1:end-1), mstxt(1:end-1), goodtxt(1:end-1));
    end
end


%*************************************************************************
function gdata = show_plot(gdata, type)

if (ischar(type))
    type = {type};
end

c = gdata.chan;
ax = gdata.axes;
d = gdata.data;

goodchan = gdata.data.goodchan(c);
med = nanmedian(d.sig(:,c));
if (any(ismember(type, {'all','plot'})))
    cla(ax,'reset');
    plot(ax, d.t,d.sig(:,c)-med,'k-', 'HitTest','off');
    axis(ax, 'tight');
    set(ax,'ButtonDownFcn',@on_click_axes);

    xl = get(ax,'XLim');    
    if (diff(xl) > 10)
        fac = diff(xl)/10;
        zoom(ax,'xon');
        zoom(ax,fac);
        zoom(ax,'off');
        zoom(ax,'on');
        zoom(ax,'off');
    end
end

if (any(ismember(type, {'all','spikes'})))
    if (ishandle(gdata.hspikes))
        delete(gdata.hspikes);
    end
    if goodchan
        gdata.hspikes = addplot(ax, d.spiket{c},d.spikeamp{c}-med, 'ro', ...
            'MarkerFaceColor','r','MarkerSize',4, 'HitTest','off');
    end
    
    if ~isempty(gdata.eventtimes)
        yl = get(ax,'YLim');
        addplot(ax, repmat(gdata.eventtimes(:)',[2 1]), ...
            repmat(yl(:),[1 length(gdata.eventtimes)]), 'y--', ...
            'HitTest','off');
    end
end

if (any(ismember(type, {'all','bursts'})))
    if (any(ishandle(gdata.hbursts)))
        delete(gdata.hbursts(ishandle(gdata.hbursts)));
    end
    
    if (goodchan)
        b = get_burst_override(gdata, c);
        
        on1 = b.on;
        off1 = b.off;

        yy = nanmean(abs(d.spikeamp{c}-med));

        h1 = addplot(ax, [on1; off1], repmat([yy; yy],[1 length(on1)]), 'b-', ...
            'LineWidth',2, 'ButtonDownFcn',@on_click_burst, ...
            'UIContextMenu', gdata.burstContextMenu);
        gdata.hbursts = double(h1);
        
        %color the overridden bursts
        set(h1(b.isover), 'Color','m');
        
        if strcmp(type,'all')
            sel = 'none';
        else
            sel = 'keep';
        end
            
        gdata = update_burst_sel(gdata,sel);

        if (ishandle(gdata.hdiag))
            gdata = show_diagnostics(gdata);
        end
    end
end

if (any(ismember(type, {'all','plot'})))
    xl = [min(d.t) max(d.t)];
    gdata.hthreshln = addplot(ax, xl,gdata.thresh(1,[c c]),'g--', ...
        xl,gdata.thresh(2,[c c]),'g--', ...
        'ButtonDownFcn',@on_click_thresh_line);
end

%*************************************************************************
function b = get_burst_override(gdata, c)

b = gdata.data.burst(c);
for i = 1:length(c)
    j = c(i);
    
    if any(~isnan(gdata.burstoverride(j).on))
        on1 = gdata.data.burst(j).on;
        off1 = gdata.data.burst(j).off;
        ctr1 = gdata.data.burst(j).ctr;
        nspike1 = gdata.data.burst(j).nspike;
        
        % add in the overrides
        on2 = gdata.burstoverride(j).on;
        off2 = gdata.burstoverride(j).off;
        
        isover1 = ~isnan(on2);
        on1(end+1:length(on2)) = NaN;
        off1(end+1:length(off2)) = NaN;
        nspike1(end+1:length(off2)) = NaN;
        
        on1(isover1) = on2(isover1);
        off1(isover1) = off2(isover1);
        ctr1(isover1) = (off2(isover1) + on2(isover1))/2;
        
        %remove bursts
        if any(gdata.burstoverride(j).remove)
            on1(gdata.burstoverride(j).remove) = NaN;
            off1(gdata.burstoverride(j).remove) = NaN;
        end

        b(i).on = on1;
        b(i).off = off1;
        b(i).ctr = ctr1;
        b(i).nspike = nspike1;
    else
        isover1 = false(size(gdata.data.burst(j).on));
    end
    b(i).isover = isover1;
end

%*************************************************************************
function gdata = update_spikes(gdata, doall)

d = gdata.data;
if (doall)
    c = 1:size(d.sig,2);
else
    c = gdata.chan;
end

for i = c
    spikeind = findspikes(d.sig(:,i), gdata.thresh(:,i));
    d.spiket{i} = d.t(spikeind{1});
    d.spikeamp{i} = d.sig(spikeind{1},i);
end
gdata.data = d;
gdata = update_bursts(gdata, doall);

%*************************************************************************
function gdata = update_bursts(gdata, doall)

d = gdata.data;
if (doall)
    c = 1:size(d.sig,2);
else
    c = gdata.chan;
end

for i = c
    if (length(d.spiket{i}) > 10)
        [burst1,spike1] = findbursts(d.spiket{i}, 'simple', 'interburstdur',gdata.interburst(i), ...
            'minspikes',gdata.minspikes(i));
        if isempty(burst1.ctr) && ~isempty(burst1.on)
            burst1.ctr = NaN;
        end
        d.burst(i) = burst1;
                
        gdata.burstoverride(i) = struct('on',NaN(1,length(burst1.on)), ...
            'off',NaN(1,length(burst1.on)), 'remove',false(1,length(burst1.on)));
    else
        d.burst(i).method = '';
        d.burst(i).interburstdur = [];
        d.burst(i).minspikes = [];
        d.burst(i).on = [];
        d.burst(i).off = [];
        d.burst(i).ctr = [];
        d.burst(i).nspike = [];
        gdata.burstoverride(c) = struct('on',[],'off',[],'remove',[]);
    end
end
gdata.data = d;

%*************************************************************************
function on_click_burst(obj, event)

gdata = guidata(obj);
s = find(gdata.hbursts == obj);
if numel(s) == 1
    gdata = update_burst_sel(gdata, s);
end
guidata(obj,gdata);

%*************************************************************************
function gdata = update_burst_sel(gdata, s)

if ischar(s)
    switch s
        case 'none'
            s = [];
        case 'keep'
            s = gdata.selburstind;
    end
end
gdata.selburstind = s;

good = ishandle(gdata.hselburst);
if any(good)
    delete(gdata.hselburst(good));
    gdata.hselburst(1:end) = -1;
end
if ~isempty(s)
    xd = get(gdata.hbursts(s),'XData');
    yd = get(gdata.hbursts(s),'YData');
    gdata.hselburst(1) = line('Parent',gdata.axes,'XData',xd,'YData',yd, ...
        'Marker','none','LineWidth',4, 'Color','m', ...
        'ButtonDownFcn',{@on_click_selected_burst,'ctr'}, ...
        'UIContextMenu', gdata.burstContextMenu);
    gdata.hselburst(2) = line('Parent',gdata.axes,'XData',xd(1),'YData',yd(1), ...
        'Marker','s','MarkerSize',10, 'Color','m', 'MarkerFaceColor','w', ...
        'ButtonDownFcn',{@on_click_selected_burst,'on'}, ...
        'UIContextMenu', gdata.burstContextMenu);
    gdata.hselburst(3) = line('Parent',gdata.axes,'XData',xd(2),'YData',yd(2), ...
        'Marker','s','MarkerSize',10, 'Color','m', 'MarkerFaceColor','w', ...
        'ButtonDownFcn',{@on_click_selected_burst,'off'}, ...
        'UIContextMenu', gdata.burstContextMenu);
    
    if ishandle(gdata.hburstrate)
        chan = gdata.chan;
        xdiag = get(gdata.hburstrate(chan),'XData');
        ydiag1 = get(gdata.hburstrate(chan),'YData');
        gdata.hselburst(4) = line('Parent',gdata.hdiag(1), 'XData', xdiag(s), ...
            'YData',ydiag1(s), 'Marker','*','MarkerSize',14, 'Color','m', ...
            'LineWidth',3);
        
        ydiag2 = get(gdata.hburstdur(chan),'YData');
        gdata.hselburst(5) = line('Parent',gdata.hdiag(2), 'XData', xdiag(s), ...
            'YData',ydiag2(s), 'Marker','*','MarkerSize',14, 'Color','m', ...
            'LineWidth',3);        
    else
        gdata.hselburst(4:5) = -1;
    end
end

%*************************************************************************
function on_click_selected_burst(obj, event, what)

gdata = guidata(obj);
c = get(gdata.axes, 'CurrentPoint');
sel = get(gdata.figure, 'SelectionType');

if strcmp(sel,'normal')
    assert(~isempty(gdata.selburstind));
    
    xd = get(gdata.hselburst(1), 'XData');
    orig = xd - c(1,1);
    set(gdata.hbursts(gdata.selburstind), 'Visible','off');
    set(gdata.figure, 'WindowButtonMotionFcn',{@on_drag_selected_burst,what,orig}, ...
        'WindowButtonUpFcn',{@on_button_up_selected_burst,what,orig});
end

%*************************************************************************
function on_drag_selected_burst(obj, event, what, orig)

gdata = guidata(obj);

c = get(gdata.axes, 'CurrentPoint');
x = c(1,1);
xd = get(gdata.hselburst(1), 'XData');

switch what
    case 'on'
        set(gdata.hselburst(1), 'XData', [x xd(2)]);
        set(gdata.hselburst(2), 'XData', x);
    case 'off'
        set(gdata.hselburst(1), 'XData', [xd(1) x]);
        set(gdata.hselburst(3), 'XData', x);

    case 'ctr'
        xd = orig + x; 
        set(gdata.hselburst(1), 'XData', xd);
        set(gdata.hselburst(2), 'XData', xd(1));
        set(gdata.hselburst(3), 'XData', xd(2));
end

%*************************************************************************
function on_button_up_selected_burst(obj, event, what,orig)

gdata = guidata(obj);

c = get(gdata.axes, 'CurrentPoint');
x = c(1,1);
xd = get(gdata.hselburst(1), 'XData');

switch what
    case 'on'
        xd = [x xd(2)];
        set(gdata.hselburst(1), 'XData', xd);
        set(gdata.hselburst(2), 'XData', x);
    case 'off'
        xd = [xd(1) x];
        set(gdata.hselburst(1), 'XData', xd);
        set(gdata.hselburst(3), 'XData', x);

    case 'ctr'
        xd = orig + x;
        set(gdata.hselburst(1), 'XData', xd);
        set(gdata.hselburst(2), 'XData', xd(1));
        set(gdata.hselburst(3), 'XData', xd(2));
end

xd = sort(xd);

s = gdata.selburstind;
set(gdata.hbursts(s), 'XData',xd, 'Color','m', 'Visible','on');

gdata.burstoverride(gdata.chan).on(s) = xd(1);
gdata.burstoverride(gdata.chan).off(s) = xd(2);

gdata = show_plot(gdata,'bursts');

set(gdata.figure, 'WindowButtonMotionFcn',[], ...
    'WindowButtonUpFcn',[]);

guidata(obj,gdata);

%*************************************************************************
function on_click_add_burst(obj,event)

gdata = guidata(obj);
c = get(gdata.axes, 'CurrentPoint');
x = c(1,1);

if ~isempty(gdata.hbursts)
    xdall = get(gdata.hbursts,'XData');
    xdall = cat(1,xdall{:})';
    d = nanmean(diff(xdall));
    
    yd = get(gdata.hbursts(1),'YData');
else
    yd = nanmean(abs(gdata.data.spikeamp{c}-med));
    yd = [yd yd];
    d = 0.1;
    xdall = [];
end

xd = x + 0.5*[-d d];

if ~isempty(xdall)
    %find the next burst after the added one
    k = first(xdall(1,:) >= x);
    
    %check for overlaps
    if isempty(k)
        if xd(1) < xdall(2,end)
            xd(1) = xdall(2,end) + 0.1*d;
        end
    elseif k == 1
        if xd(2) > xdall(1,1)
            xd(2) = xdall(1,1) - 0.1*d;
        end
    elseif ~isempty(k) && (k > 1)
        if xd(1) < xdall(2,k-1)
            xd(1) = xdall(2,k-1) + 0.1*d;
        end
        if xd(2) > xdall(1,k)
            xd(2) = xdall(1,k) - 0.1*d;
        end
    end
end

gdata.burstoverride(gdata.chan).on(end+1) = xd(1);
gdata.burstoverride(gdata.chan).off(end+1) = xd(2);
gdata.burstoverride(gdata.chan).remove(end+1) = false;

hb = line('Parent',gdata.axes, 'XData',xd, 'YData',yd, 'Color','m', ...
    'LineWidth',2, 'ButtonDownFcn',@on_click_burst, ...
    'UIContextMenu', gdata.burstContextMenu);
gdata.hbursts(end+1) = hb;

if (ishandle(gdata.hdiag))
    gdata = show_diagnostics(gdata);
end

guidata(obj,gdata);


%*************************************************************************
function on_click_split_burst(obj,event)

gdata = guidata(obj);
c = get(gdata.axes, 'CurrentPoint');

xdall = get(gdata.hbursts,'XData');
xdall = cat(1,xdall{:})';

ind = find((c(1,1) >= xdall(1,:)) & (c(1,1) <= xdall(2,:)), 1);
if ~isempty(ind)
    xd0 = xdall(:,ind);

    gap = diff(xd0)/10;
    dur = (diff(xd0)-gap)/2;
    
    xd1 = xd0(1) + [0 dur];
    xd2 = xd0(2) + [-dur 0];
    
    gdata.burstoverride(gdata.chan).on(end+1) = xd2(1);
    gdata.burstoverride(gdata.chan).off(end+1) = xd2(2);
    gdata.burstoverride(gdata.chan).remove(end+1) = false;
    
    h1 = gdata.hbursts(ind);
    h2 = copyobj(h1, gdata.axes);
    set(h1,'XData',xd1);
    set(h2,'XData',xd2, 'ButtonDownFcn',@on_click_burst, ...
        'UIContextMenu', gdata.burstContextMenu);
    
    gdata.hbursts(end+1) = h2;
    
    if (ishandle(gdata.hdiag))
        gdata = show_diagnostics(gdata);
    end
    gdata = update_burst_sel(gdata,ind);
    
    guidata(obj,gdata);
end




%*************************************************************************
function on_click_remove_burst(obj,event)

gdata = guidata(obj);
c = get(gdata.axes, 'CurrentPoint');

xdall = get(gdata.hbursts,'XData');
xdall = cat(1,xdall{:})';

ind = find((c(1,1) >= xdall(1,:)) & (c(1,1) <= xdall(2,:)), 1);
if ~isempty(ind)
    gdata.burstoverride(gdata.chan).on(ind) = 0;
    gdata.burstoverride(gdata.chan).off(ind) = 0;
    gdata.burstoverride(gdata.chan).remove(ind) = true;

    set(gdata.hbursts(ind),'XData',[NaN; NaN], 'YData',[NaN; NaN]);

    if (ishandle(gdata.hdiag))
        gdata = show_diagnostics(gdata);
    end
    gdata = update_burst_sel(gdata,'none');

    guidata(obj,gdata);
end


%*************************************************************************
function on_click_thresh_line(obj, event)

gdata = guidata(obj);
yd = get(obj,'YData');
set(gdata.figure, 'WindowButtonMotionFcn',@(o,e) on_drag_thresh_line(o,e,sign(yd(1))), ...
    'WindowButtonUpFcn',@(o,e) on_button_up_thresh_line(o,e,sign(yd(1))));

%*************************************************************************
function on_drag_thresh_line(obj, event, s)

gdata = guidata(obj);

c = get(gdata.axes, 'CurrentPoint');
y = c(1,2);
if (s > 0)
    if (y > 0)
        set(gdata.hthreshln(2), 'YData',[y y]);
    end
else
    if (y < 0)
        set(gdata.hthreshln(1), 'YData',[y y]);
    end
end


%*************************************************************************
function on_button_up_thresh_line(obj, event, s)

gdata = guidata(obj);

c = get(gdata.axes, 'CurrentPoint');
y = c(1,2);
good = false;
if ((s > 0) && (y > 0))
    set(gdata.hthreshln(2), 'YData',[y y]);
    gdata.thresh(2,gdata.chan) = y;
    good = true;
elseif ((s < 0) && (y < 0))
    set(gdata.hthreshln(1), 'YData',[y y]);
    gdata.thresh(1,gdata.chan) = y;
    good = true;
end

if (good)
    if (s < 0)
        set(gdata.spikeThreshLoEdit, 'String', num2str(y,3));
    else
        set(gdata.spikeThreshHiEdit, 'String', num2str(y,3));
    end
    
    gdata = update_spikes(gdata, false);
    gdata = show_plot(gdata, {'spikes','bursts'});
end

set(gdata.figure, 'WindowButtonMotionFcn',[], ...
    'WindowButtonUpFcn',[]);

guidata(obj,gdata);

%*************************************************************************
function on_set_chan(obj,event,dchan,chan)

gdata = guidata(obj);
chan0 = gdata.chan;

if (~isempty(chan))
    gdata.chan = chan;
else
    gdata.chan = gdata.chan + dchan;
end

nchan = size(gdata.data.sig,2);
if ((gdata.chan < 1) || (gdata.chan > nchan))
    gdata.chan = chan0;
else
    if (gdata.chan == nchan)
        set(gdata.nextButton,'Enable','off');
        set(gdata.prevButton,'Enable','on');
    elseif (gdata.chan == 1)
        set(gdata.nextButton,'Enable','on');
        set(gdata.prevButton,'Enable','off');
    else
        set(gdata.nextButton,'Enable','on');
        set(gdata.prevButton,'Enable','on');
    end
    
    gdata = show_plot(gdata,'all');
end
set(gdata.channelEdit, 'String', num2str(gdata.chan));
set(gdata.spikeThreshLoEdit, 'String', num2str(gdata.thresh(1,gdata.chan),3));
set(gdata.spikeThreshHiEdit, 'String', num2str(gdata.thresh(2,gdata.chan),3));
set(gdata.interburstDurEdit, 'String', num2str(gdata.interburst(gdata.chan),3));
set(gdata.minSpikesEdit, 'String', num2str(gdata.minspikes(gdata.chan),3));
set(gdata.skipChannelCheck, 'Value', ~gdata.data.goodchan(gdata.chan));

guidata(obj,gdata);

%*************************************************************************
function on_edit_channel(obj,event)

gdata = guidata(obj);
s = get(obj,'String');
c = str2double(s);
guidata(obj,gdata);

on_set_chan(obj,event,[],c);

%*************************************************************************
function on_edit_spikeThresh(obj,event,sgn)

gdata = guidata(obj);
s = get(obj,'String');
c = str2double(s);

good = false;
if (sgn < 0) && (c < 0)
    gdata.thresh(1,gdata.chan) = c;
    good = true;
elseif (sgn > 0) && (c > 0)
    gdata.thresh(2,gdata.chan) = c;
    good = true;
end
if (good)
    gdata = update_spikes(gdata, false);
    gdata = show_plot(gdata, {'spikes','bursts'});
else
    if (sgn < 0)
        set(obj,'String',num2str(gdata.thresh(1,gdata.chan)));
    else
        set(obj,'String',num2str(gdata.thresh(2,gdata.chan)));
    end        
end
guidata(obj,gdata);


%*************************************************************************
function on_edit_interburst(obj,event)

gdata = guidata(obj);
s = get(obj,'String');
c = str2double(s);

gdata.interburst(gdata.chan) = c;
gdata = update_bursts(gdata, false);
gdata = show_plot(gdata, 'bursts');

guidata(obj,gdata);


%*************************************************************************
function on_edit_minSpikes(obj,event)

gdata = guidata(obj);
s = get(obj,'String');
c = str2double(s);

gdata.minspikes(gdata.chan) = c;
gdata = update_bursts(gdata, false);
gdata = show_plot(gdata, 'bursts');

guidata(obj,gdata);

%*************************************************************************
function on_click_burstdiagnostics(obj,event)

gdata = guidata(obj);
if (get(obj, 'Value') > 0)
    pos = get(gdata.figure,'Position');
    pos(1) = pos(1) + pos(3);
    f = figure('WindowStyle','normal', 'Tag','burstrate');
    set(f,'Position',pos);
    
    gdata.hdiag(1) = subplot(2,1,1,'Parent',f);
    ylabel(gdata.hdiag(1), 'Burst rate (Hz)');
    gdata.hdiag(2) = subplot(2,1,2,'Parent',f);
    ylabel(gdata.hdiag(2), 'Burst duration (sec)');
    xlabel(gdata.hdiag(2), 'Time (sec)');
    linkaxes(gdata.hdiag,'x');
    
    marker = 'osdp^v';
    col = 'bgrmc';
    nchan = size(gdata.data.sig,2);
    
    if length(marker) < nchan
        r = ceil(nchan/length(marker));
        marker = repmat(marker,[1 r]);
    end
    if length(col) < nchan
        r = ceil(nchan/length(col));
        col = repmat(col,[1 r]);
    end
    gdata.diagmarker = marker;
    gdata.diagcol = col;
    
    gdata.hburstrate = zeros(nchan,1);
    gdata.hburstdur = zeros(nchan,1);
    for i = 1:nchan
        gdata.hburstrate(i) = line('XData',[],'YData',[], 'Color',col(i), ...
            'Marker',marker(i), 'LineStyle','none', 'Parent',gdata.hdiag(1), ...
            'ButtonDownFcn',{@on_click_burst_diag, gdata.figure});
        gdata.hburstdur(i) = line('XData',[],'YData',[], 'Color',col(i), ...
            'Marker',marker(i), 'LineStyle','none', 'Parent',gdata.hdiag(2), ...
            'ButtonDownFcn',{@on_click_burst_diag, gdata.figure});
    end
    
    gdata = show_diagnostics(gdata);
else
    f = findobj('Tag','burstrate');
    delete(f);
    
    gdata.hburstrate = -1;
    gdata.hburstdur = -1;
end
guidata(obj,gdata);

%*************************************************************************
function gdata = show_diagnostics(gdata)

chan = find(gdata.data.goodchan);

b = get_burst_override(gdata,chan);

for i = 1:length(b)
    ctr1 = b(i).ctr;
    [ctr1,ord] = sort(ctr1);
    burstfreq1 = NaN(size(ctr1));
    burstfreq1(ord(2:end-1)) = 2./(ctr1(3:end) - ctr1(1:end-2));
    
    j = chan(i);
    gdata.data.burstfreq{j} = burstfreq1;
    gdata.data.burstdur{j} = b(i).off - b(i).on;
end

for i = 1:length(chan)
    j = chan(i);
    
    m = gdata.diagmarker(j);
    c = gdata.diagcol(j);
    s = [];
    if j == gdata.chan
        c = 'k';
        fc = 'k';
        sz = 12;
        if ~isempty(gdata.selburstind)
            s = gdata.selburstind;
        end
    else
        fc = 'none';
        sz = 8;
    end

    set(gdata.hburstrate(j), 'XData',b(i).ctr, ...
        'YData',gdata.data.burstfreq{j}, 'Marker',m, 'Color',c, ...
        'MarkerFaceColor',fc, 'MarkerSize',sz);
    set(gdata.hburstdur(j), 'XData',b(i).ctr, ...
        'YData',gdata.data.burstdur{j}, 'Marker',m, 'Color',c, ...
        'MarkerFaceColor',fc, 'MarkerSize',sz);
    
    if ~isempty(s) && ishandle(gdata.hselburst(4)) && ishandle(gdata.hselburst(5))
        set(gdata.hselburst(4), 'XData', b(i).ctr(s), ...
            'YData',gdata.data.burstfreq{j}(s));
        set(gdata.hselburst(5), 'XData', b(i).ctr(s), ...
            'YData',gdata.data.burstdur{j}(s));
    end    
end
axis(gdata.hdiag(1), 'tight');
axis(gdata.hdiag(2), 'tight');

if ~isempty(gdata.eventtimes)
    yl = get(gdata.hdiag(1),'YLim');
    addplot(gdata.hdiag(1), repmat(gdata.eventtimes(:)',[2 1]), ...
        repmat(yl(:),[1 length(gdata.eventtimes)]), 'y--', ...
        'HitTest','off');
    
    yl = get(gdata.hdiag(2),'YLim');
    addplot(gdata.hdiag(2), repmat(gdata.eventtimes(:)',[2 1]), ...
        repmat(yl(:),[1 length(gdata.eventtimes)]), 'y--', ...
        'HitTest','off');
end


%*************************************************************************
function on_click_burst_diag(obj,event, fig)

gdata = guidata(fig);

ax = get(obj,'Parent');
c = get(ax, 'CurrentPoint');
chan = find((gdata.hburstrate == obj) | (gdata.hburstdur == obj));
if numel(chan) == 1
    if chan ~= gdata.chan
        on_set_chan(fig,event,[],chan);
        gdata = guidata(fig);
    end

    ctrx = c(1,1);
    ctry = c(1,2);
    xl = get(gdata.axes, 'XLim');
    dx = diff(xl)/2;
    set(gdata.axes, 'XLim', ctrx + [-dx dx]);
    
    xdata = get(obj, 'XData');
    ydata = get(obj, 'YData');

    xl2 = get(ax, 'XLim');
    dx2 = diff(xl2);
    yl2 = get(ax, 'YLim');
    dy2 = diff(yl2);
    
    [~,ind] = min(((ctrx - xdata)/dx2).^2 + ((ctry - ydata)/dy2).^2);
    gdata = update_burst_sel(gdata, ind);
end

guidata(fig,gdata);

    
%*************************************************************************
function on_click_skipchannel(obj,event)

gdata = guidata(obj);
if (~get(obj,'Value'))
    gdata.data.goodchan(gdata.chan) = true;
    gdata = update_spikes(gdata, false);
else
    gdata.data.goodchan(gdata.chan) = false;
    i = gdata.chan;
    gdata.data.spiket{i} = [];
    gdata.data.spikeamp{i} = [];
    gdata = update_bursts(gdata, false);
end  

show_plot(gdata,{'busts','spikes'});
guidata(obj,gdata);

%*************************************************************************
function on_click_axes(obj,event)

gdata = guidata(obj);

sel = get(gdata.figure, 'SelectionType');

switch sel
    case 'normal'
        % turn off the selected burst
        gdata = update_burst_sel(gdata, 'none');
    case 'open'
        on_click_add_burst(obj,event);
        gdata = guidata(obj);
end
guidata(obj,gdata);


%*************************************************************************
function on_click_Done(obj,event)

data = guidata(obj);
uiresume(data.figure);

%*************************************************************************
function setUICallbacks(data)

set(data.doneButton,'Callback',@on_click_Done);
set(data.nextButton, 'Callback',{@on_set_chan,1,[]});
set(data.prevButton, 'Callback',{@on_set_chan,-1,[]});
set(data.channelEdit, 'Callback',@on_edit_channel);
set(data.spikeThreshLoEdit, 'Callback',{@on_edit_spikeThresh,-1});
set(data.spikeThreshHiEdit, 'Callback',{@on_edit_spikeThresh,1});
set(data.interburstDurEdit, 'Callback',@on_edit_interburst);
set(data.minSpikesEdit, 'Callback',@on_edit_minSpikes);
set(data.burstDiagnosticsCheck, 'Callback',@on_click_burstdiagnostics);
set(data.skipChannelCheck, 'Callback',@on_click_skipchannel);

set(data.axes, 'ButtonDownFcn',@on_click_axes);

set(data.axes, 'UIContextMenu', data.axesContextMenu);
set(data.addBurstMenu, 'Callback', @on_click_add_burst);

set(data.splitBurstMenu, 'Callback', @on_click_split_burst);
set(data.removeBurstMenu, 'Callback', @on_click_remove_burst);











