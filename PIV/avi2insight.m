function avi2insight(aviname,outdir,varargin)

opt.camera = 'L';
opt.traverse = 0;
opt.delay = 0;
opt.process = 0;
opt.hardwareconfig = 0;
opt.pairskip = 1;
opt.fieldskip = 1;
opt.frames = [1 Inf];
opt.basename = '';
opt.filespec = '%s%06d.T%03d.D%03d.P%03d.H%03d.%c%c.tif';

opt.display = true;

opt = parsevarargin(opt,varargin, 3);

if nargin == 0
    [fn,pn] = uigetfile({'*.avi';'*.cine';'*.*'},'Select movie file');
    aviname = fullfile(pn,fn);
    vr = VideoReader(aviname);
    nframes = vr.NumberOfFrames;
    
    outdir = uigetdir(pn, 'Select output directory');
    
    [~,basename,~] = fileparts(fn);
    
    prompt = {'Base file name:',...
        'Pair skip:',...
        'Field skip:',...
        'Start frame:',...
        'End frame:',...
        'Frame rate:',...
        'Display:'};
    def = {basename, '1','1','1',num2str(nframes),'300','on'};
    
    vals = inputdlg(prompt,'Options',1,def);
    
    basename = vals{1};
    opt.pairskip = str2double(vals{2});
    opt.fieldskip = str2double(vals{3});
    opt.frames = [str2double(vals{4}), str2double(vals{5})];
    opt.framerate = str2double(vals{6});
    if ismember(lower(vals{7}),{'on','yes','y'})
        opt.display = true;
    else
        opt.display = false;
    end;
else
    vr = VideoReader(aviname);
    nframes = vr.NumberOfFrames;
    if isempty(opt.basename)
        [~,basename,~] = fileparts(aviname);
    else
        basename = opt.basename;
    end
end

if opt.frames(2) > nframes
    opt.frames(2) = nframes;
end

aframes = opt.frames(1):opt.fieldskip:opt.frames(2);
bframes = aframes + opt.pairskip;

good = (aframes >= 1) & (bframes <= nframes);
aframes = aframes(good);
bframes = bframes(good);
nconv = sum(good);

if opt.display
    fig = figure;
    aax = axes('Position',[0 0 0.5 1]);
    bax = axes('Position',[0.5 0 0.5 1]);

    afr = read(vr, aframes(1));
    aim = image(afr(:,:,1), 'Parent',aax);
    axis(aax, 'tight','equal','ij','off');
    bfr = read(vr, bframes(1));
    bim = image(bfr(:,:,1), 'Parent',bax);
    axis(bax, 'tight','equal','ij', 'off');
    
    colormap gray;
else
    timedWaitBar(0, 'Saving frames...');
end

for i = 1:nconv
    afr = read(vr, aframes(i));
    aname = sprintf(opt.filespec, basename, aframes(i), opt.traverse, opt.delay, ...
        opt.process, opt.hardwareconfig, opt.camera,'A');
    
    bfr = read(vr, bframes(i));
    bname = sprintf(opt.filespec, basename, aframes(i), opt.traverse, opt.delay, ...
        opt.process, opt.hardwareconfig, opt.camera,'B');
    
    if opt.display
        set(aim,'CData',afr(:,:,1));
        set(bim,'CData',bfr(:,:,1));
        drawnow;
    end
    
    imwrite(afr(:,:,1), fullfile(outdir,aname), 'tif');
    imwrite(bfr(:,:,1), fullfile(outdir,bname), 'tif');
    
    if ~opt.display
        timedWaitBar(i/nconv);
    end
end

close(fig);


        
    


    