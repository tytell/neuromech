function [t, imnum, dv, raw] = getPCOtime(vid, varargin)

opt.frames = [];
opt = parsevarargin(opt,varargin,2);

[pn,fn,ext] = fileparts(vid);
if (strcmpi(ext,'.tif'))
    tifinfo = imfinfo(vid);
    
    nfr = numel(tifinfo);
    if (isempty(opt.frames))
        frames = 1:nfr;
    else
        frames = opt.frames;
    end
    
    timedWaitBar(0,'Getting timestamps');
    
    imnum = zeros(nfr,1);
    dv = zeros(nfr,6);
    raw = zeros(nfr,14,'uint8');
    for i = frames
        I = imread(vid, i, 'Info',tifinfo);
        
        [dv(i,:),imnum(i),raw(i,:)] = getPCOtimestamp(I);
        timedWaitBar(i/nfr);
    end
    timedWaitBar(1);
    
    t = dv(:,4)*3600 + dv(:,5)*60 + dv(:,6);
end

        
        
    