function h5imagestack(filenames,outname,varargin)

opt.deflate = 5;
opt.chunk = [16 16 16];
opt.crop = [];
opt.scale = [];     % micron/voxel
opt.checkcrop = true;

opt = parsevarargin(opt,varargin,3);

I = imread(filenames{1});
imshow(I);
addplot(opt.crop([1 2 2 1 1]),opt.crop([3 3 4 4 3]), 'r-');
if ~inputyn('Is crop OK? ','default',true)
    disp('Canceling');
    return;
end

%transpose so that x refers to the first dimension and y to the second
I = I';

if ~isempty(opt.crop)
    I = I(opt.crop(1):opt.crop(2),opt.crop(3):opt.crop(4));
end

sz = size(I);

N = length(filenames);

dt = class(I);
h5create(outname,'/Image',[sz N], 'Datatype',dt, 'ChunkSize',opt.chunk, ...
    'Deflate',opt.deflate);

rng = [Inf -Inf];

timedWaitBar(0,'Loading images');
for i = 1:N
    I = imread(filenames{i});
    I = I';
    if ~isempty(opt.crop)
        I = I(opt.crop(1):opt.crop(2),opt.crop(3):opt.crop(4));
    end
    h5write(outname,'/Image',I, [ones(size(sz)) i], [sz 1]);
    
    rng(1) = min([rng(1); I(:)]);
    rng(2) = max([rng(2); I(:)]);
    timedWaitBar(i/N);
end
timedWaitBar(1);

if ~isempty(opt.scale)
    h5writeatt(outname,'/Image','Scale',opt.scale);
end
h5writeatt(outname,'/Image','Minimum',rng(1));
h5writeatt(outname,'/Image','Maximum',rng(2));





