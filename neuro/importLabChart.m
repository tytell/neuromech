function importLabChart(filename, outname, varargin)

opt.channelnames = 'auto';
opt = parsevarargin(opt, varargin);

F = load(filename);

if isfield(F,'datastart')
    %new file type
    nchan = length(F.datastart);
    ismultirate = any(F.samplerate ~= F.samplerate(1));
    if (ismultirate)
        t = cell(1,nchan);
    end
    data = cell(1,nchan);
    for i = 1:nchan
        data{i} = F.data(F.datastart(i):F.dataend(i))';

        if (ismultirate)
            t{i} = (0:length(data{i})-1)'/F.samplerate(i) + F.firstsampleoffset(i);
        end
    end
    if (~ismultirate)
        t{1} = (0:length(data{1})-1)'/F.samplerate(1);
    end
else
    assert(isfield(F,'data_block1'))
    nchan = size(F.data_block1,1);
    data = mat2cell(F.data_block1', size(F.data_block1,2), ones(1,nchan));
    t{1} = F.ticktimes_block1';
    ismultirate = false;
end

if isfield(F,'comtext_block1')
    comtext = mat2cell(F.comtext_block1,ones(size(F.comtext_block1,1),1),size(F.comtext_block1,2));
else
    comtext = mat2cell(F.comtext,ones(size(F.comtext,1),1),size(F.comtext,2));
end
evttxt = 'Event Marker';
evtind = find(strncmp(comtext, evttxt, length(evttxt)), 1, 'first');

if isfield(F,'com')
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
    iscom0(evtind) = false;

    iscom = iscom0(F.com(:,5));
    commenttxt = comtext(F.com(iscom,5));
    commentt = F.com(iscom,3) / F.tickrate;

    isspike = isspike0(F.com(:,5));
    spiket = F.com(isspike,3) / F.tickrate;
    spikenum = spikenum0(F.com(isspike,5));
    spikeunit = spikeunit0(F.com(isspike,5));
    spikeheight = spikeheight0(F.com(isspike,5));
    spikewidth = spikewidth0(F.com(isspike,5));
    isspike = true;
else
    isspike = false;
end

if ischar(opt.channelnames)
    if isfield(F,'titles_block1')
        ttl = 'titles_block1';
    else
        ttl = 'titles';
    end
    switch opt.channelnames
        case 'auto'
            channelnames = cell(1,nchan);
            for i = 1:nchan
                if goodchan(i)
                    channelnames{i} = genvarname(F.titles(i,:));
                end
            end
    end
elseif iscell(opt.channelnames)
    channelnames = opt.channelnames;
end

S = cell2struct(data(goodchan),channelnames(goodchan),2);

if ismultirate
    for i = 1:nchan
        if goodchan(i)
            tname1 = genvarname([channelnames{i} 't']);
            S.(tname1) = t{i};
        end
    end
else
    S.t = t{1};
end

if isspike
    S.spiket = spiket;
    S.spikenum = spikenum;
    S.spikeheight = spikeheight;
    S.spikewidth = spikewidth;
    S.spikeunit = spikeunit;
    S.eventt = eventt;
end
S.eventt = eventt;
S.blocktimes = F.blocktimes;
S.tickrate = F.tickrate;
S.samplerate = F.samplerate;
S.commentt = commentt;
S.commenttxt = commenttxt;
if isfield(F,'tickrate')
    S.blocktimes = F.blocktimes;
    S.tickrate = F.tickrate;
    S.samplerate = F.samplerate;
    S.commentt = commentt;
    S.commenttxt = commenttxt;
end

save(outname,'-struct','S');



