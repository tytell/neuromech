function [t, imnum, dv] = getPCOtime(vid, varargin)

opt.frames = [];
opt.quiet = false;
opt = parsevarargin(opt,varargin,2);

[~,~,ext] = fileparts(vid);
if (strcmpi(ext,'.tif'))
    tifinfo = imfinfo(vid);
    
    nfr = numel(tifinfo);
    if (isempty(opt.frames))
        frames = 1:nfr;
    else
        frames = opt.frames;
    end
    
    if (~opt.quiet)
        timedWaitBar(0,'Getting timestamps');
    end
    
    imnum = zeros(nfr,1);
    dv = zeros(nfr,6);
    for i = frames
        I = imread(vid, i, 'Info',tifinfo);
        
        [dv(i,:),imnum(i)] = getPCOtimestamp(I);
        if (~opt.quiet)
            timedWaitBar(i/nfr);
        end
    end
    if (~opt.quiet)
        timedWaitBar(1);
    end;
    
    t = dv(:,4)*3600 + dv(:,5)*60 + dv(:,6);
end

        
        
    