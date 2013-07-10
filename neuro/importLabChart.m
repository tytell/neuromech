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
        t{i} = (0:length(data{i})-1)*F.samplerate(i) + F.firstsampleoffset(i);
    end
end
if (~ismultirate)
    t{1} = (0:length(data{1})-1)*F.samplerate(1);
end

comtext = mat2cell(F.comtext,ones(size(comtext,1)),size(comtext,2));
tok = regexp(comtext, 'Spike: height = ([\d.]+) V, width = ([\d.]+) ms, Spike (\d+), Unit (\d+)', ...
    'tokens','once');
    

