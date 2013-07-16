function importLabChart(filename, outname, varargin)

opt.channelnames = 'auto';
opt = parsevarargin(opt, varargin);

F = load(filename);

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

comtext = mat2cell(F.comtext,ones(size(F.comtext,1),1),size(F.comtext,2));
tok = regexp(comtext, 'Spike: height = ([\d.]+) V, width = ([\d.]+) ms, Spike (\d+)(, Unit )?(\d+)?', ...
    'tokens','once');
isspike = cellfun(@(x) ~isempty(x), tok);

spiket = t{1}(F.com(isspike,3));
spikenum = cellfun(@(x) str2double(x{3}), tok(isspike));
spikeheight = cellfun(@(x) str2double(x{1}), tok(isspike));
spikewidth = cellfun(@(x) str2double(x{2}), tok(isspike));
spikeunit = cellfun(@(x) str2double(x{5}), tok(isspike));

if ischar(opt.channelnames)
    switch opt.channelnames
        case 'auto'
            channelnames = cell(1,nchan);
            for i = 1:nchan
                channelnames{i} = genvarname(F.titles(i,:));
            end
    end
elseif iscell(opt.channelnames)
    channelnames = opt.channelnames;
end

S = cell2struct(data,channelnames,2);

if ismultirate
    for i = 1:nchan
        tname1 = genvarname([channelnames{i} 't']);
        S.(tname1) = t{i};
    end
else
    S.t = t{1};
end

S.spiket = spiket;
S.spikenum = spikenum;
S.spikeheight = spikeheight;
S.spikewidth = spikewidth;
S.spikeunit = spikeunit;

save(outname,'-struct','S');



