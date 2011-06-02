function varargout = findspikes(varargin)
% function [ind,thresh] = findspikes(y, opts...)
%  or                     findspikes(t,y, opts...)
%
% Locates spikes in a time series, based on a threshold.  Optionally asks
% for the threshold

opt.nspike = 2;
opt.window = false;
opt.method = 'threshold';
opt.threshold = [];
opt.maxspikedur = 0.02;     % in sec
opt.mininterspike = 0.001;  % sec
opt.nsmooth = 15;
opt.subtractdc = true;
opt.highpasscutoff = 500;   % Hz
opt.debug = false;
opt.downsampledisplay = 1;

if ((nargin == 1) || ((nargin >= 2) && ischar(varargin{2})))
    p = 2;
else
    p = 3;
end;

opt = parsevarargin(opt,varargin(p:end),p, 'multival',{...
    'method',{'threshold','intracell'}});

%sort out the threshold if necessary
switch opt.method,
    case 'threshold',
        if (p == 2)
            sig = varargin{1};
            threshval = opt.threshold;
        elseif (p == 3)
            [sig,threshval] = varargin{1:2};
            if (isempty(threshval))
                threshval = opt.threshold;
            end;
        end;
        
        sig = shiftdim(sig);
        nchan = size(sig,2);
        
        if (isempty(threshval)),
            %assume 5000Hz and take the first 10 sec
            if (size(sig,1) < 100000),
                k = 1:size(sig,1);
            else
                k = 1:100000;
            end;
            
            box = msgbox(['Click on axes to set threshold for each channel.' ...
                ' Hit any key to end.'],...
                'Threshold','help','non-modal');
            
            %run through each channel and define the threshold
            for i = 1:nchan,
                clf;
                zoom off;
                
                h0 = plot(k,sig(k,i)-nanmedian(sig(k,i)));
                
                set(h0,'HitTest','off');
                thresh0 = max(abs(prctile(sig(i,:),[15 85])));
                
                h = addplot([k(1) k(1); k(end) k(end)],...
                    [thresh0 -thresh0; thresh0 -thresh0],'r-');
                set(h,'HitTest','off');
                set(gca,'ButtonDownFcn',{@setThresh,h});
                set(gcf,'KeyPressFcn',@endThresh, 'UserData','');
                
                waitfor(gcf,'UserData');
                
                thresh0 = get(h(1),'YData');
                threshval(i) = thresh0(1);
                delete(h);
            end;
            if (ishandle(box)),
                delete(box);
            end;
            set(gca,'ButtonDownFcn','');
            set(gcf,'KeyPressFcn','');
        end;
        
        ind = cell(nchan,1);
        for chan = 1:nchan,
            %subtract the median and rectify
            asig = abs(sig(:,chan) - nanmedian(sig(:,chan)));
            len = length(asig);
            
            %find spikes
            %first, find anything above the spike threshold
            spike1 = asig >= threshval(chan);
            %now make sure they're real peaks, by checking that each one is the
            %maximum value in the nspike elements around it
            k = opt.nspike+1:len-opt.nspike;
            for i = 1:opt.nspike,
                spike1(k) = spike1(k) & (asig(k+i) < asig(k)) & ...
                    (asig(k-i) < asig(k));
            end;
            ind{chan} = find(spike1);
        end;
        
        if (nargout == nchan+1),
            %NB: this means that for a single channel, the index will be returned
            %as a vector, *not* as a cell
            varargout = ind;
            varargout{nchan+1} = threshval;
        elseif (nargout == nchan),
            varargout = ind;
        elseif (nargout == 2),
            varargout = {ind, threshval};
        else
            varargout = {ind};
        end;

    case 'intracell',
        if (p == 2)
            error('findspikes:params','Must pass both a time and a signal for intracell spike finding');
        end;
        
        [dt,sig] = varargin{1:2};
        if (numel(dt) > 1)
            dt = dt(2) - dt(1);
        end;
        
        sig = makecol(sig);        
        nchan = size(sig,2);
        
        mininterspike = opt.mininterspike/dt;

        if (opt.subtractdc)
            [b,a] = butter(5,opt.highpasscutoff*dt/2, 'low');
        
            siglow = filtfilt(b,a,sig);
        
            sig = sig - siglow;
        end;
        
        if (isempty(opt.threshold) || isnan(opt.threshold))
            opt.threshold = thresholddlg(dt,sig,'max', 'downsample',opt.downsampledisplay);
        end;
        
        [pk,pkind] = findpeaks2(sig, 'threshold',opt.threshold, ...
            'minpeakdistance', mininterspike, 'bias','biggerpeaks');
        
        if (opt.debug),
            t = (0:length(sig)-1)*dt;
            
            plot(t,sig, t(pkind),pk,'r.');
        end;
        
        switch nargout
            case 1,
                varargout = {pkind};
            case 2,
                varargout = {pkind, opt.threshold};
        end;
end;
                
            
        

%------------------------------------------------------------------------
function setThresh(obj,eventdata, hLine)

pt = get(obj,'CurrentPoint');

set(hLine(1),'YData',abs([pt(1,2) pt(1,2)]));
set(hLine(2),'YData',-abs([pt(1,2) pt(1,2)]));


%------------------------------------------------------------------------
function endThresh(obj,eventdata)

set(obj,'UserData','done');
