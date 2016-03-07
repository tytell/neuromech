function varargout = importLabChart(filename, outname, varargin)

opt.channelnames = 'auto';
opt.outformat = 'old';
opt = parsevarargin(opt, varargin);

switch opt.outformat
    case 'old'
        isoldformat = true;
    case 'new'
        isoldformat = false;
    otherwise
        error('Unrecognized output format');
end

F = load(filename);

if ~isfield(F,'datastart')
    % "Simple" format
    data = num2cell(F.data_block1,2);
    data = data';
    nchan = length(data);
    for i = 1:nchan
        data{i} = data{i}';
    end
    t = {F.ticktimes_block1'};
    ismultirate = false;
    
    commenttxt1 = mat2cell(F.comtext_block1,ones(size(F.comtext_block1,1),1),size(F.comtext_block1,2));
    commentind1 = F.comtick_block1;
    commentt1 = t{1}(F.comtick_block1);
    commentchan1 = F.comchan_block1;
    
    commentt = cell(nchan,1);
    commentind = cell(nchan,1);
    commenttxt = cell(nchan,1);
    for i = 1:nchan
        ischan = commentchan1 == i;
        commentt{i} = commentt1(ischan);
        commentind{i} = commentind1(ischan);
        commenttxt{i} = commenttxt1(ischan);
    end

    titles = F.titles_block1;
    goodchan = 1:nchan;
else
    goodchan = find(F.datastart >= 1);
    nchan = length(goodchan);
    ismultirate = any(F.samplerate(goodchan) ~= F.samplerate(goodchan(1)));
    if (ismultirate)
        t = cell(1,nchan);
    end
    data = cell(1,nchan);
    for i = 1:nchan
        j = goodchan(i);
        data{i} = F.data(F.datastart(j):F.dataend(j))';

        if (ismultirate)
            t{i} = (0:length(data{i})-1)'/F.samplerate(j) + F.firstsampleoffset(j);
        end
    end
    if (~ismultirate)
        t{1} = (0:length(data{1})-1)'/F.samplerate(goodchan(1));
    end

    comtext = mat2cell(F.comtext,ones(size(F.comtext,1),1),size(F.comtext,2));
    evttxt = 'Event Marker';
    evtind = find(strncmp(comtext, evttxt, length(evttxt)), 1, 'first');

    if isempty(evtind)
        evtind = -1;
    end
    isevt = F.com(:,5) == evtind;
    eventt = F.com(isevt,3) / F.tickrate;

    tok = regexp(comtext, 'Spike: height = ([\d.]+) V, width = (-?[\d.]+) ms, Spike (\d+)(, Unit )?(\d+)?', ...
        'tokens','once');

    spikenum0 = NaN(size(comtext));
    spikeheight0 = NaN(size(comtext));
    spikewidth0 = NaN(size(comtext));
    spikeunit0 = NaN(size(comtext));
    for i = 1:length(tok)
        if (~isempty(tok{i}))
            spikenum0(i) = str2double(tok{i}{3});
            spikeheight0(i) = str2double(tok{i}{1});
            spikewidth0(i) = str2double(tok{i}{2});
            spikeunit0(i) = str2double(tok{i}{5});
        end
    end
    isspike0 = ~isnan(spikenum0);

    iscom0 = ~isspike0;
    if evtind > 0
        iscom0(evtind) = false;
    end

    iscom = iscom0(F.com(:,5));
    commenttxt = comtext(F.com(iscom,5));
    commentt = F.com(iscom,3) / F.tickrate;

    isspike = isspike0(F.com(:,5));
    spiket = F.com(isspike,3) / F.tickrate;
    spikenum = spikenum0(F.com(isspike,5));
    spikeunit = spikeunit0(F.com(isspike,5));
    spikeheight = spikeheight0(F.com(isspike,5));
    spikewidth = spikewidth0(F.com(isspike,5));
    
    titles = F.titles;
end

if ischar(opt.channelnames)
    switch opt.channelnames
        case 'auto'
            channelnames = cell(1,nchan);
            for i = 1:nchan
                channelnames{i} = genvarname(titles(goodchan(i),:));
            end
    end
elseif iscell(opt.channelnames)
    channelnames = opt.channelnames;
end

if isoldformat
    S = cell2struct(data,channelnames,2);
else
    S.channels = cell2struct(data,channelnames,2);
end

if ismultirate
    for i = 1:nchan
        if isoldformat
            tname1 = genvarname([channelnames{i} 't']);
            S.(tname1) = t{i};
        else
            S.t.(channelnames{i}) = t{i};
        end
    end
else
    S.t = t{1};
end

if (exist('spiket','var'))
    S.spiket = spiket;
    S.spikenum = spikenum;
    S.spikeheight = spikeheight;
    S.spikewidth = spikewidth;
    S.spikeunit = spikeunit;
    S.eventt = eventt;
    S.blocktimes = F.blocktimes;
end
if (isfield(F,'tickrate'))
    S.tickrate = F.tickrate;
    S.samplerate = F.samplerate;
end
S.commentt = commentt;
S.commentind = commentind;
S.commenttxt = commenttxt;

if ~isempty(outname)
    save(outname,'-struct','S');
else
    varargout = {S};
end



