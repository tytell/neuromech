function threshold = getthreshold(t,sig,varargin)
% function threshold = getthreshold(t,sig,options...)

threshold = [];
nthresh = 1;
channames = {};
showdur = 10;

p = 1;
while (p <= length(varargin)),
    switch lower(varargin{p}),
        case 'nthresh',
            nthresh = varargin{p+1};
            p = p+2;
            
        case 'threshold',
            threshold = varargin{p+1};
            p = p+2;

        case 'channelnames',
            channames = varargin{p+1};
            p = p+2;
        
        case 'showdur',
            showdur = varargin{p+1};
            p = p+2;
            
        otherwise,
            error('Unrecognized option %s\n',varargin{p});
    end;
end;

%get the signal data structured right
t = makecol(t);
if (size(sig,2) > size(sig,1)),
    sig = sig';
end;
nchan = size(sig,2);

%bring the current figure to the front
figure(gcf);

if (~isempty(threshold)),
    nthresh = size(threshold,1);
end;

if (nthresh == 1),
    box = msgbox(['Click on axes to set threshold for each channel.' ...
        ' Hit any key to end.'],...
        'Threshold','help','non-modal');
elseif (nthresh == 2),
    box = msgbox(['Left (right) click on axes to set high (low) threshold for each channel.' ...
        ' Hit any key to end.'],...
        'Threshold','help','non-modal');
end;

if (isempty(threshold)),
    threshold = nans(nthresh,nchan);
end;

ispos = all(sig(:) >= 0);
isshow = t <= showdur;
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
    thresh0 = threshold(:,i);
    if ((nthresh == 1) && isnan(thresh0)),
        thresh0 = max(abs(prctile(sig(:,i),[15 85])));
    elseif (nthresh == 2),
        thresh00 = [max(abs(prctile(sig(:,i),[5 95]))) max(abs(prctile(sig(:,i),[15 85])))];
        undef = isnan(thresh0);
        thresh0(undef) = thresh00(undef);
    end;

    %plot the threshold line(s)
    if (~ispos),
        x = repmat([t(1); maxt], [1 2*length(thresh0)]);
        y = repmat(thresh0', [4 1]);
        y = reshape(y,[2 2*length(thresh0)]);
        if (nthresh == 2),
            y(:,[2 4]) = -y(:,[2 4]);
        else
            y(:,2) = -y(:,2);
        end;
    else
        x = repmat([t(1); maxt], [1 length(thresh0)]);
        y = repmat(thresh0', [2 1]);
        y = reshape(y,[2 length(thresh0)]);
    end;
    h = addplot(x,y,'r-');
    set(h(3:end),'Color','g');
    
    set(h,'HitTest','off');
    set(gca,'ButtonDownFcn',{@localSetThresh,h});
    set(gcf,'KeyPressFcn',@localEndThresh, 'UserData','');
    if (~isempty(channames)),
        title(sprintf('Channel %d (%s)',i,channames{i}));
    end;

    waitfor(gcf,'UserData');

    if (nthresh == 1),
        th = get(h(1),'YData');
        thresh0 = th(1);
    else
        th = get(h(1),'YData');
        thresh0(1) = th(1);
        th = get(h(3),'YData');
        thresh0(2) = th(1);
    end;
    threshold(:,i) = thresh0';
    delete(h);
end;
if (ishandle(box)),
    delete(box);
end;
set(gca,'ButtonDownFcn','');
set(gcf,'KeyPressFcn','');


%------------------------------------------------------------------------
function localSetThresh(obj,eventdata, hLine)

pt = get(obj,'CurrentPoint');
fig = get(obj,'Parent');
sel = get(fig,'SelectionType');

if (length(hLine) == 1),
    i = 1;
    ispos = true;
elseif (length(hLine) == 2),
    i = 1;
    ispos = false;
elseif (length(hLine) == 4),
    switch sel,
        case 'normal',
            i = 1;
        case 'alt',
            i = 3;
        otherwise,
            i = 1;
    end;
    ispos = false;
end;

set(hLine(i),'YData',abs([pt(1,2) pt(1,2)]));
if (~ispos),
    set(hLine(i+1),'YData',-abs([pt(1,2) pt(1,2)]));
end;

%------------------------------------------------------------------------
function localEndThresh(obj,eventdata)

set(obj,'UserData','done');
