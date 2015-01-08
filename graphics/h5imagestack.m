function h5imagestack(filenames,outname,varargin)

opt.deflate = 5;
opt.chunk = [16 16 16];
opt.crop = [];
opt.scale = [];     % micron/voxel

opt = parsevarargin(opt,varargin,3);

I = imread(filenames{1});
if ~isempty(opt.crop)
    I = I(opt.crop(3):opt.crop(4), opt.crop(1):opt.crop(2));
end

sz = size(I);

N = length(filenames);


dt = class(I);
h5create(outname,'/Image',[sz N], 'Datatype',dt, 'ChunkSize',opt.chunk, ...
    'Deflate',opt.deflate);

timedWaitBar(0,'Loading images');
for i = 1:N
    I = imread(filenames{i});
    if ~isempty(opt.crop)
        I = I(opt.crop(3):opt.crop(4), opt.crop(1):opt.crop(2));
    end
    h5write(outname,'/Image',I, [ones(size(sz)) i], [sz 1]);
    timedWaitBar(i/N);
end
timedWaitBar(1);

if ~isempty(opt.scale)
    h5writeatt(outname,'/Image','Scale',opt.scale);
end


