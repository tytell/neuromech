function varargout = printNeuralTrace(varargin)
%function printNeuralTrace(t,sig,...)
%   or    printNeuralTrace(filename.daq,...)
%
% Options: 'spikes', [style] - Find spikes (if necessary) and plot them without the
%              raw traces.  Can pass an additional spike style: 'raw'
%              (default), which just plots whatever it gets; 'up', plots
%              spikes only going up; 'upconst', plots constant height
%              spikes going up; 'updown', plots symmetric spikes going up
%              and down; 'updownconst', same as updown, but with constant
%              height; 'bar', plots spikes without a zero line (like a bar
%              code); 'barconst', like bar but constant height
%          'spikethreshold' - Threshold for finding spikes
%          'bursts', burstpos - Plots burst position defined by burst pos,
%              which can be either a two column set of spike indices
%              (defining the beginning and end of each burst) or a column
%              vector of burst center times
%          'margins', [left bottom width height] - Page margins in inches
%          'lineheight', height - Line height in inches (default 0.75")
%          'linegap', gap - Gap between lines in inches (default 0.08")
%          'linedur', dur - Time duration to shown in a line (default
%              10sec)
%          'channeloffset', offset - Offset between channels within a line
%              in relative height (the default, 1, means no overlap between
%              channels.  0.75 would overlap by 25%, and so forth)
%          'channelnames', names - Names of the channels.  Not necessary if
%              the data are read from a .daq file.
%          'stimchan', chan - Index of the channel with the stimulus trace
%              (default 1)
%          'nostim' - Do not plot a stimulus trace
%          'alignstimulus' - Arrange successive lines so that stimulus
%              periods line up down the page
%          'stimbars' - Show gray and white bars to indicate positive or
%              negative stimulus times
%          'col' or 'color' - Sets colors for the channels.  Can be a
%              single color character, or enough for all the channels
%          'noprint' - Do not automatically print the page
%          'append' - Append the trace to an existing printNeuralTrace plot


%check whether they passed in a .daq file name or just the trace itself
if ((isnumeric(varargin{1}) || iscell(varargin{1})) && ...
        ((nargin == 1) || ~isa(varargin{2},class(varargin{1})))),
    t = varargin{1};
    if (iscell(t)),
        for i = 1:length(t),
            sig{i} = ones(size(t{i}));
        end;
    end;
    p = 2;
    isFile = false;
elseif ((nargin >= 2) && (isnumeric(varargin{1}) || iscell(varargin{1})) && ...
        isa(varargin{2},class(varargin{1}))),
    t = varargin{1};
    sig = varargin{2};
    p = 3;
    isFile = false;
elseif ((nargin >= 1) && ischar(varargin{1}) && ~isempty(varargin{1})),
    filename = varargin{1};
    p = 2;
    isFile = true;
elseif ((nargin == 0) || isempty(varargin{1})),
    [fn,pn] = uigetfile('*.daq','Choose data acquisition file');
    filename = fullfile(pn,fn);
    p = 2;
    isFile = true;
end;
opts = varargin(p:end);

%default options
findSpikes = false;
isSpikes = false;
showBursts = false;
spikethreshold = [];
margins = [0.5 0.5 0.75 1];         % left right top bottom in inches
lineheight = 0.75;                  % height of a line of traces in inches
linegap = 0.08;                     % gap between lines in inches
linedur = 10;                       % duration of a single line in sec
channeloffset = 1;                  % offset between channels in relative height
channelnames = {};
channels = [];
stimchan = 1;
isStimulus = true;
alignStimulus = false;              % align successive traces to stimulus period
stimBars = false;
spikeStyle = 'raw';
spikeStyles = {'raw','up','upconst','updown','updownconst','bar','barconst'};
isAxis = true;
cols = 'bgrcyk';
printopts = {};
isPrint = true;
isAppend = false;

%process all the options
i = 1;
while (i <= length(opts)),
    switch lower(opts{i}),
        case {'spike','spikes'},
            findSpikes = true;
            if ((i < length(opts)) && ischar(opts{i+1}) && ...
                 ~isempty(strmatch(lower(opts{i+1}),spikeStyles,'exact'))),
                spikeStyle = opts{i+1};
                i = i+2;
            else
                i = i+1;
            end;
            
        case {'spikethreshold','spikethresh'},
            spikethreshold = opts{i+1};
            i = i+2;
            
        case 'bursts',
            burstval = opts{i+1};
            showBursts = true;
            i = i+2;
            
        case 'margins',
            margins = opts{i+1};
            i = i+2;
            
        case 'lineheight',
            lineheight = opts{i+1};
            i = i+2;
            
        case 'linegap',
            linegap = opts{i+1};
            i = i+2;
            
        case 'linedur',
            linedur = opts{i+1};
            i = i+2;
            
        case {'channeloffset','chanoff'},
            channeloffset = opts{i+1};
            i = i+2;
            
        case {'channel','channels'},
            channels = opts{i+1};
            i = i+2;
            
        case 'channelnames',
            channelnames = opts{i+1};
            i = i+2;
            
        case {'stimchan','stimuluschannel'},
            stimchan = opts{i+1};
            i = i+2;
            
        case {'nostim','nostimulus'},
            isStimulus = false;
            i = i+1;
            
        case 'alignstimulus',
            alignStimulus = true;
            i = i+1;
            
        case 'stimbars',
            stimBars = true;
            i = i+1;
            
        case 'noAxes',
            isAxis = false;
            i = i+1;

        case {'col','color','colors'},
            cols = opts{i+1};
            i = i+2;
            
        case 'noprint',
            isPrint = false;
            i = i+1;
            
        case 'append',
            isAppend = true;
            i = i+1;
            
        otherwise,
            %assume its an option for the print function
            printopts{end+1} = opts{i};
            i = i+1;
    end;
end;

%load in the daq file if that's what they passed in
if (isFile),
    info = daqread(filename,'info');
    nSamples = info.ObjInfo.SamplesAcquired;
    channelnames = {info.ObjInfo.Channel.ChannelName};
    daqdate = datestr(info.ObjInfo.InitialTriggerTime, 0);

    %get the stimulus
    if (~isempty(stimchan)),
        stimname = channelnames{stimchan};
        [stimt,stimphase,stim] = loadStimulus(filename, 'Channels',stimchan, ...
            'Downsample',100);
        stim = (stim - nanmedian(stim))/range(stim);
    else
        stimname = '';
        isStimulus = false;
    end;
    if (~isempty(channels)),
        channels = channels(channels <= length(channelnames));
        channelnames = channelnames(channels);
    else
        channels = 1:length(info.ObjInfo.Channel);
    end;
    nChannels = length(channels);

    if (findSpikes),
        %find the spikes
        if (isempty(spikethreshold)),
            spikefig = figure;
        else
            spikefig = [];
        end;
        [sig,t,spikethreshold] = loadSpikes(filename, 'Channels',channels, ...
            'Threshold',spikethreshold);
        if (~isempty(spikefig)),
            delete(spikefig);
        end;
        isAbsSig = false;
        isSpikes = true;
    else
        %or read the raw data
        [sig,t] = daqread(filename, 'Channels',channels);
        isAbsSig = false;
    end;
else
    %process the raw trace
    if (iscell(sig)),
        nSamples = max(cellfun(@length,sig));
        nChannels = length(sig);     
    else
        if (size(sig,1) < size(sig,2)),
            warning('printNeuralTrace:dataOrientation',...
                'Data appears to be arranged along rows.  Transposing...');
            sig = sig';
            if (size(t,1) > 1),
                t = t';
            end;
        end;
        
        nSamples = size(sig,1); 
        nChannels = size(sig,2);
    end;
    
    %first pull the stimulus out
    if (isStimulus),
        nChannels = nChannels - 1;
        
        if (iscell(sig)),
            rest = [1:stimchan-1 stimchan+1:length(sig)];
            stim = shiftdim(sig{stimchan});
            sig = sig(rest);
            stimt = shiftdim(t{stimchan});
            t = t(rest);
        else
            rest = [1:stimchan-1 stimchan+1:size(sig,2)];
            stim = sig(:,stimchan);
            sig = sig(:,rest);
            stimt = t(:,stimchan);
            if (size(t,2) > 1),
                t = t(:,rest);
            end;
        end;
        if (~isempty(channelnames)),
            stimname = channelnames(stimchan);
            channelnames = channelnames(rest);
        else
            stimname = 'stimulus';
        end;
        stimphase = getStimPhase(stimt,stim);
    end;
    
    %check for irregular spacing in time - that means they passed spikes
    if (iscell(t)),
        t1 = t{1};
        sig1 = sig{1};
    else
        t1 = t(:,1);
        sig1 = sig(:,1);
    end;
    dt1 = t1(2)-t1(1);
    if (any(abs(diff(t1) - dt1)/dt1 > 0.01)),
        isSpikes = true;
    end;
    if (all(sig1 >= 0)),
        isAbsSig = true;
    else
        isAbsSig = false;
    end;
    
    %set up channel names if we don't have them already
    if (isempty(channelnames)),
        for i = 1:nChannels,
            channelnames{i} = sprintf('ch%d',i);
        end;
    end;
    filename = '';
end;

%normalize the channels
chanheight = zeros(nChannels,1);
for i = 1:nChannels,
    if (iscell(sig)),
        sig1 = shiftdim(sig{i});
    else
        sig1 = sig(:,i);
    end;
    
    chanheight(i) = 2*prctile(abs(sig1),95);
    chanheight(i) = min(max(abs(sig1)),chanheight(i));
    
    if (~isAbsSig),
        sig1 = sig1 - median(sig1);
    end;
    sig1 = sig1/chanheight(i);
    
    if (iscell(sig)),
        sig{i} = sig1;
        
        t{i} = shiftdim(t{i});
    else
        sig(:,i) = sig1;
    end;
end;

%check what sort of burst data we have
if (~showBursts),
    burstType = 'none';
elseif ((size(burstval,2) == 2) && all(burstval(:) == floor(burstval(:)))),
    %two columns and they're all integers -> indices
    burstind = burstval;
    burstType = 'index';
elseif (any(size(burstval) == 1)),
    %vector and not integers -> time values
    bursttime = burstval;
    burstType = 'time';
end;

%deal with the colors
if (isStimulus),
    stimcol = cols(end);
    if (length(cols) > 1),
        cols = cols(1:end-1);
    end;
    
    stim = stim - nanmedian(stim);
    stim = stim/range(stim) * channeloffset;
else
    alignStimulus = false;
end;
if (length(cols) < nChannels),
    cols = repmat(cols,[1 ceil(nChannels/length(cols))]);
end;

%display the figure on screen at 1.5x smaller than the page printout
figReduce = 1.5;

%create the figure and set up the scaling correctly
if (isPrint || isAppend),
    fig = findobj('Tag','printNeuralTraceFigure');
else
    fig = [];
end;
if (isempty(fig)),
    fig = figure('Tag','printNeuralTraceFigure');
else
	fig = fig(1);
    figure(fig);
end;
figs = fig;

if (isAppend),
    %if we're appending then get the existing axes
    hAx = findobj(fig, 'Tag','printNeuralTraceAxis');
    hHead = findobj(fig, 'Tag','printNeuralTraceHeader');
    hFoot = findobj(fig, 'Tag','printNeuralTraceFooter');
    
    if (isempty(hAx)),
        isAppend = false;
    else
        nAxes = length(hAx);
        %sort the axes by position
        pos = get(hAx,'Position');
        pos = cat(1,pos{:});
        %sort by increasing bottom position
        [q,ord] = sort(pos(:,2));
        %and reverse to make sure we have the top axis first
        hAx = hAx(ord(end:-1:1));

        visible = strmatch('on',get(hAx, 'Visible'));
        appendFirstAxis = visible(end)+2;
        if (appendFirstAxis > length(hAx)),
            isAppend = false;
        end;
    end;
end;

%sometimes they want us to append, but we can't, so double check whether
%isAppend changed above
if (~isAppend),
    clf;
    
    %set up the figure to have the same proportions as the page size
    set(fig,'PaperUnits','inches','Units','inches');
    papersz = get(fig,'PaperSize');
    paperpos = [0 0 papersz];
    figpos = get(fig,'Position');
    %keep the top of the figure in the same place
    figpos(2) = figpos(2) + figpos(4) - papersz(2)/figReduce;
    %but update it to have a half page size with the right aspect ratio
    figpos(3:4) = papersz/figReduce;

    %and set it to fill the page
    set(fig, 'Position',figpos, 'PaperPositionMode','manual',...
        'PaperPosition',paperpos, 'Color','w');

    %figure out how many axes are going to fit
    printsz = papersz - [margins(1)+margins(2) margins(3)+margins(4)];
    nAxes = floor(printsz(2) / (lineheight + linegap));
    linegap = (printsz(2) - (nAxes * lineheight))/(nAxes - 1);

    %create all the axes.
    top = papersz(2) - margins(3);
    hAx = zeros(nAxes,1);
    for i = 1:nAxes,
        axpos = [margins(1) top - i*lineheight - (i-1)*linegap ...
            printsz(1) lineheight];
        axpos = axpos/figReduce;
        hAx(i) = axes('Units','inches', 'Color','none', 'LineWidth',1, ...
                      'Tag','printNeuralTraceAxis');
        set(hAx(i), 'Position',axpos);
        set(hAx(i), 'Units','normalized');
    end;

    %create header and footer axes
    hHead = axes('Units','inches', 'Color','none', 'LineWidth',1, ...
                 'Tag','printNeuralTraceHeader');
    axpos = [margins(1) papersz(2)-margins(3) printsz(1) margins(3)/2];
    axpos = axpos/figReduce;
    set(hHead, 'Position',axpos);
    set(hHead, 'Units','normalized','Visible','off');

    hFoot = axes('Units','inches', 'Color','none', 'LineWidth',1,...
                 'Tag','printNeuralTraceFooter');
    axpos = [margins(1) margins(4)/2 printsz(1) margins(4)/2];
    axpos = axpos/figReduce;
    set(hFoot, 'Position',axpos);
    set(hFoot, 'Units','normalized','Visible','off');
    
    appendFirstAxis = 1;
end;

%get the range of times
if (iscell(t)),
    tinit = min(cellfun(@min,t));
    tfinal = max(cellfun(@max,t));
else
    tinit = min(t(:));
    tfinal = max(t(:));
end;
if (isStimulus),
    tinit = min(tinit,min(stimt));
    tfinal = max(tfinal,max(stimt));
end;

%and figure out how many total lines we'll need, then how many pages those
%lines will fit on
nLines = ceil((tfinal-tinit)/linedur);
nPages = ceil((nLines - (nAxes-appendFirstAxis+1))/nAxes) + 1;

%now run through the pages
tcur = tinit;
linedur(1) = linedur;
for page = 1:nPages,
    %on the first page, we'll do the appending if necessary
    if (page == 1),
        ax1 = appendFirstAxis;
    else
        ax1 = 1;
    end;
    lastAxis = nAxes+1;
    
    %now run through all the axes on this page
    for ax = ax1:nAxes,
        axes(hAx(ax));
        cla(hAx(ax),'reset');
        set(hAx(ax),  'Color','none', 'LineWidth',1, ...
                      'Tag','printNeuralTraceAxis');
                  
        if (tcur <= tfinal),
            set(hAx(ax),'Visible','on');
        else
            %make extra axes invisible
            lastAxis = min(lastAxis,ax-1);
            set(hAx(ax),'Visible','off');
            continue;
        end;

        hold on;

        %if we're aligning the stimulus, figure out how much we need to
        %overlap with the previous line to make them line up
        if (alignStimulus),
            %first find the times for the current line
            k = find((stimt >= tcur) & (stimt < tcur+linedur));
            %then look for points when the stimulus phase resets back to
            %zero
            stimperind = k(diff(stimphase(k)) < 0);
            if (length(stimperind) > 1),
                %use those to estimate the stimulus period
                stimper = median(diff(stimt(stimperind)));

                %we should overlap with the previous line by however much
                %of a period extra we'll need to make the first period
                %start on an integer multiple of the stimulus period
                alignOverlap = stimper - rem(stimt(stimperind(1))-tcur, stimper);
                tcur = tcur - alignOverlap;
            end;
        end;

        %make sure we're showing the whole duration
        lim = [tcur tcur+linedur 0 channeloffset*nChannels];
        %set the vertical limits appropriately
        if (isStimulus),
            lim(3) = lim(3) - channeloffset;
        elseif (~isAbsSig && ...
                isempty(strmatch(lower(spikeStyle),{'up','upconst'},'exact'))),
            lim(3:4) = lim(3:4) - channeloffset/2;
        end;

        %now plot the bars for the stimulus
        if (isStimulus && stimBars),
            k = find((stimt >= tcur) & (stimt < tcur+linedur));
            
            if (any(abs(stim(k)) > 0.1)),
                %points where the phase is between 0 and pi get a bar - but to
                %make the calculations smoother (i.e., to avoid having to deal
                %with the wrap in the phase value at 2 pi), we add on pi/2
                phaseplus = mod(stimphase(k) + pi/2, 2*pi);
                %look for points where the phase passes through pi/2 or starts
                %up from a NaN
                barleft = k((isnan(phaseplus(1:end-2)) | (phaseplus(1:end-2) < pi/2)) & ...
                    (phaseplus(2:end-1) >= pi/2) & isfinite(phaseplus(3:end)));
                %same in reverse
                barright = k((phaseplus(1:end-2) < 3*pi/2) & ...
                    (phaseplus(2:end-1) >= 3*pi/2) & isfinite(phaseplus(3:end)));

                if (~isempty(barleft) && ~isempty(barright)),
                    %check to make sure the first left edge is before the first
                    %right edge
                    if (barright(1) < barleft(1)),
                        a = first(isfinite(stimphase(k)));
                        barleft = [k(a); barleft];
                    end;
                    %and the last right edge should be greater than the last
                    %left edge
                    if (barleft(end) > barright(end)),
                        b = last(isfinite(stimphase(k)));
                        barright = [barright; k(b)];
                    end;

                    %if barleft and barright are different lengths, then
                    %something's screwed up - attempt to recover
                    if (length(barleft) ~= length(barright)),
                        q = 1:min(length(barleft),length(barright));
                        barleft = barleft(q);
                        barright = barright(q);
                    end;

                    %plot the bars
                    barx = [stimt(barleft) stimt(barleft) stimt(barright) stimt(barright)]';
                    bary = repmat([lim([3 4]) lim([4 3])]',[1 length(barleft)]);
                    fill(barx,bary,[0.8 0.8 0.8], 'EdgeColor','none');
                end;
            end;
        end;

        %plot the channel names in the first axis
        if ((ax == ax1) && ~isempty(channelnames{1})),
            minleft = tcur+linedur;
            for chan = 1:nChannels,
                off = channeloffset*nChannels - chan*channeloffset;
                if (~isAbsSig),
                    off = off + channeloffset/2;
                end;
                h = text(tcur+linedur,off, channelnames{chan}, ...
                    'HorizontalAlignment','right','VerticalAlignment','middle', ...
                    'Color',cols(chan));
                ext = get(h,'Extent');
                
                %careful because the text font size has to be the print
                %size, even though the figure is reduced.  So we have to
                %divide the extent by figReduce to get the actual print
                %extent
                minleft = min(minleft,tcur+linedur - ext(3)/figReduce);
            end;
            if (isStimulus && ~isempty(stimname)),
                off = off - 0.5*channeloffset;
                h = text(tcur+linedur,off, stimname, ...
                    'HorizontalAlignment','right','VerticalAlignment','middle',...
                    'Color',stimcol);
                ext = get(h,'Extent');
                minleft = min(minleft,tcur+linedur - ext(3)/figReduce);
            end;
            linedur1 = minleft - tcur;
        else
            linedur1 = linedur;
        end;

        %run through the channels
        for chan = 1:nChannels,
            if (iscell(t)),
                t1 = t{chan};
                sig1 = sig{chan};
            else
                if (size(t,2) == nChannels),
                    t1 = t(:,chan);
                end;
                sig1 = sig(:,chan);
            end;
            k = find((t1 >= tcur) & (t1 < tcur+linedur1));

            off = channeloffset*nChannels - chan*channeloffset;
            if (~isAbsSig),
                off = off + channeloffset/2;
            end;

            %plot spikes
            if (isSpikes),
                switch lower(spikeStyle),
                    case 'raw',
                        tt = repmat(t1(k)',[3 1]);
                        ss = [zeros(size(k)) sig1(k) zeros(size(k))]';
                        
                    case 'up',
                        tt = repmat(t1(k)',[3 1]);
                        ss = [zeros(size(k)) abs(sig1(k)) zeros(size(k))]';
                        
                    case 'upconst',
                        tt = repmat(t1(k)',[3 1]);
                        ss = [zeros(size(k)) ones(size(k)) zeros(size(k))]';
                        
                    case 'updown',
                        tt = repmat(t1(k)',[4 1]);
                        ss = [zeros(size(k)) abs(sig1(k)/2) -abs(sig1(k)/2) zeros(size(k))]';
                        off = off + channeloffset/2;
                        
                    case 'updownconst',
                        tt = repmat(t1(k)',[4 1]);
                        ss = [zeros(size(k)) ones(size(k)) -ones(size(k)) zeros(size(k))]';
                        off = off + channeloffset/2;
                        
                    case 'bar',
                        tt = [t1(k) t1(k) repmat(NaN,size(k))]';
                        ss = [abs(sig1(k)/2) -abs(sig1(k)/2) repmat(NaN,size(k))]';
                        off = off + channeloffset/2;
                        
                    case 'barconst',
                        tt = [t1(k) t1(k) repmat(NaN,size(k))]';
                        ss = [ones(size(k)) -ones(size(k)) repmat(NaN,size(k))]';
                        off = off + channeloffset/2;
                end;
                
                switch burstType,
                    case 'index',
                        burstind1 = burstind{chan};

                        inrng = any((burstind1 >= k(1)) & (burstind1 <= k(end)),2);
                        burstind1 = burstind1(inrng,:);
                        if (burstind1(1,1) < k(1)),
                            burstind1(1,1) = k(1);
                        end;
                        if (burstind1(end,2) > k(end)),
                            burstind1(end,2) = k(end);
                        end;

                        isburst = false(size(k));
                        for i = 1:size(burstind1,1),
                            isburst((burstind1(i,1):burstind1(i,2))-k(1)+1) = true;
                        end;

                        plot(flatten(tt(:,~isburst)),flatten(ss(:,~isburst))+off, ...
                            'Color',cols(chan), 'LineWidth',0.5, 'Clipping','off');
                        plot(tt(:,isburst),ss(:,isburst)+off, ...
                            'Color',cols(chan), 'LineWidth',2);
                        plot(t1(burstind1),zeros(size(burstind1))+off,'.', 'Color',cols(chan));
                        
                    case 'time',
                        burstt1 = bursttime{chan};
                        
                        inrng = (burstt1 >= t1(k(1))) & (burstt1 <= t1(k(end)));
                        plot(tt(:),ss(:)+off, 'Color',cols(chan), 'Clipping','off');
                        plot(burstt1(inrng),zeros(sum(inrng),1)+off,'o', ...
                            'MarkerEdgeColor','w','MarkerFaceColor',cols(chan));
                        
                    otherwise,
                        plot(tt(:),ss(:)+off, 'Color',cols(chan), 'Clipping','off');
                end;
            else
                plot(t1(k), sig1(k) + off, 'Color',cols(chan), 'Clipping','off');
            end;
        end;
        
        if (isStimulus),
            off = off - 0.5*channeloffset;
            k = find((stimt >= tcur) & (stimt < tcur+linedur1));
            plot(stimt(k), stim(k) + off, stimcol);
        end;

        axis(lim);
        if (isAxis),
            xt = get(hAx(ax), 'XTick');
            if ((ax == nAxes) || (tcur + linedur > tfinal)),
                %on the bottom axis, keep all the tick labels
                set(hAx(ax), 'YColor','w', 'YTick',[], 'XTick',xt, ...
                    'XTickMode','manual','XTickLabelMode','manual');
            else
                %otherwise only show the beginning and end
                xtlab = cell(1,length(xt));
                xtlab{1} = num2str(xt(1));
                xtlab{end} = num2str(xt(end));
                set(hAx(ax), 'YColor','w', 'YTick',[], 'XTick',xt, ...
                    'XTickMode','manual','XTickLabelMode','manual',...
                    'XTickLabel',xtlab);
            end;
        else
            axis off;
            text(tcur,(lim(3)+lim(4))/2, sprintf('%.1fsec',tcur), ...
                'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
                'Rotation',90);
        end;
        hold off;
        
        tcur = tcur + linedur1;
    end;

    %reverse the layering order of the axes so that the tick labels on
    %upper axes will be visible
    for ax = nAxes:-1:1,
        axes(hAx(ax));
    end;
    
    %make the new header and footer
    if (isAppend && (page == 1)),
        axes(hAx(appendFirstAxis-1));
        delete(allchild(gca));
    else
        axes(hHead);
        delete(allchild(hHead));
    end;

    if (isFile),
        axis([0 1 0 1]);
        [pn,fn,ext] = fileparts(filename);
        text(0,0, sprintf('%s (%d samples)',[fn ext], nSamples), ...
            'HorizontalAlignment','left','VerticalAlignment','bottom');
        text(1,0,daqdate, ...
            'HorizontalAlignment','right','VerticalAlignment','bottom');
    end;
    
    if (~isAppend && (nPages > 1)),
        axes(hFoot);
        delete(allchild(hFoot));
        text(0.5,0, sprintf('Page %d of %d',page,nPages), ...
            'HorizontalAlignment','center','VerticalAlignment','top');
    end;
    
    if (isPrint),
        if (~isempty(printopts)),
            if ((length(printopts) == 1) && any(printopts{1} == ' ')),
                po = printopts{1};
                spaces = find(po == ' ');
                spaces = [0 spaces length(po)+1];
                printopts = cell(1,length(spaces)-1);
                for i = 1:length(spaces)-1,
                    printopts{i} = po(spaces(i)+1:spaces(i+1)-1);
                end;
            end;

            print(printopts{:});
        else
            printdlg;
        end;
    elseif (page < nPages),
        figs(page) = copyobj(fig,0);
    end;
end;
figs(nPages) = fig;

out = {figs,lastAxis < nAxes-1,spikethreshold};
varargout = out(1:nargout);
