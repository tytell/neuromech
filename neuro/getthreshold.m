function threshval = getthreshold(t,sig,varargin)
% function threshval = getthreshold(t,sig,options...)

opt.threshold = [];
opt.symmetric = true;
opt.channames = {};
opt.showdur = 10;
opt = parsevarargin(opt,varargin,3);

threshval = opt.threshold;

%get the signal data structured right
t = makecol(t);
if (size(sig,2) > size(sig,1)),
    sig = sig';
end;
nchan = size(sig,2);

%bring the current figure to the front
figure(gcf);

if (size(threshval,1) == 1)
    threshval = repmat(threshval,[2 1]);
end
if (size(threshval,2) == 1)
    threshval = repmat(threshval,[1 nchan]);
end

box = msgbox(['Click on axes to set threshold for each channel.' ...
    ' Hit any key to end.'],...
    'Threshold','help','non-modal');

if (isempty(threshval)),
    threshval = NaN(2,nchan);
end;

ispos = all(sig(:) >= 0);
isshow = t <= opt.showdur;
maxt = max(t(isshow));

%run through each channel and define the threshold(s)
for i = 1:nchan,
    clf;
    zoom off;

    if (~ispos)
        sig(:,i) = sig(:,i) - nanmedian(sig(:,i));
    end;
    
    plot(t(isshow),sig(isshow,i),'HitTest','off');

    %set up the default values
    thresh0 = threshval(:,i);
    if (isnan(thresh0)),
        if opt.symmetric
            thresh0 = [-1 1]*max(abs(prctile(sig(:,i),[15 85])));
        else
            thresh0 = prctile(sig(:,i),[15 85]);
        end
    end

    %plot the threshold line(s)
    if (~ispos),
        x = repmat([t(1); maxt], [1 length(thresh0)]);
        y = repmat(thresh0, [2 1]);
    else
        x = repmat([t(1); maxt], [1 length(thresh0)]);
        y = repmat(thresh0', [2 1]);
        y = reshape(y,[2 length(thresh0)]);
    end;
    h = addplot(x,y,'r-');
    
    set(h,'HitTest','off');
    set(gca,'ButtonDownFcn',{@localSetThresh,h,opt.symmetric});
    set(gcf,'KeyPressFcn',@localEndThresh, 'UserData','');
    if (~isempty(opt.channames)),
        title(sprintf('Channel %d (%s)',i,opt.channames{i}));
    end;

    waitfor(gcf,'UserData');

    if (opt.symmetric),
        th = get(h(2),'YData');
        thresh0 = [-1;1]*th(1);
    else
        th = get(h(1),'YData');
        thresh0(1) = th(1);
        th = get(h(2),'YData');
        thresh0(2) = th(1);
    end;
    threshval(:,i) = thresh0';
    delete(h);
end;
if (ishandle(box)),
    delete(box);
end;
set(gca,'ButtonDownFcn','');
set(gcf,'KeyPressFcn','');

if opt.symmetric
    threshval = threshval(2,:);
end


%------------------------------------------------------------------------
function localSetThresh(obj,~, hLine,issym)

pt = get(obj,'CurrentPoint');
fig = get(obj,'Parent');

y1 = get(hLine(1),'YData');
y2 = get(hLine(2),'YData');

if (issym)
    yval = [-1;1] * abs(pt(1,2));
else
    if (abs(pt(1,2) - y1(1)) < abs(pt(1,2) - y2(1)))
        yval(1) = pt(1,2);
        yval(2) = y2(1);
    else
        yval(1) = y1(1);
        yval(2) = pt(1,2);
    end
end

set(hLine(1),'YData',yval([1 1]));
set(hLine(2),'YData',yval([2 2]));

%------------------------------------------------------------------------
function localEndThresh(obj,eventdata)

set(obj,'UserData','done');
