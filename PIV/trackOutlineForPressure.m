function trackOutlineForPressure

[vecfilenames,pathname] = uigetfile('*.vec','Choose vector files', ...
    'MultiSelect','on');
if isempty(pathname)
    return;
end

vectok = regexp(vecfilenames, '(.*)(\d{6}).T\d+.D\d+.P(\d+).H\d+.L.vec', 'once','tokens');

vecnum = zeros(size(vecfilenames));
Pnum = NaN(size(vecfilenames));
good = true(size(vecfilenames));

for i = 1:length(vecfilenames)    
    if length(vectok{i}) == 3
        vecnum(i) = str2double(vectok{i}{2});
        Pnum(i) = str2double(vectok{i}{3});
    else
        good(i) = false;
    end
end

vecnum = vecnum(good);
Pnum = Pnum(good);
vecfilenames = vecfilenames(good);

Pvals = unique(Pnum(~isnan(Pnum)));
if length(Pvals) > 1
    fprintf('Different P values found: ');
    fprintf('%d ', Pvals);
    fprintf('\n');
    
    useP = [];
    good = false;
    for i = 1:10
        useP = input('Use P value: ','s');
        if ~isempty(useP) && any(Pvals == str2double(useP))
            useP = str2double(useP);
            good = true;
            break;
        end
    end
    
    if ~good
        fprint('No P value chosen. Cancelling');
        return;
    end
    
    isP = Pnum == useP;
    
    vecnum = vecnum(isP);
    Pnum = Pnum(isP);
    vecfilenames = vecfilenames(isP);
end
              
basename = vectok{first(good)}{1};

[basepath,vecdir] = fileparts(pathname(1:end-1));

if ~exist(fullfile(basepath,'RawData'),'dir')
    imgdir = uigetdir(basepath, 'Find image file directory');
else
    imgdir = fullfile(basepath,'RawData');
end

imgfilenames = cell(size(vecfilenames));
good = true(size(vecfilenames));
for i = 1:length(vecfilenames)
    imgfilenames{i} = fullfile(imgdir, ...
        sprintf('%s%06d.T000.D000.P000.H000.LA.tif', basename, vecnum(i)));
    
    if ~exist(imgfilenames{i},'file')
        warning('Image file %s does not exist. Skipping', imgfilenames{i});
        good(i) = false;
    end
end

if any(~good)
    if ~inputyn(sprintf('%d image files skipped (%d%%). Continue?', sum(~good), ...
            round(sum(~good) / length(good) * 100)), 'default', true)
        fprintf('Cancelling...');
        return;
    end
end

vecfilenames = vecfilenames(good);
imgfilenames = imgfilenames(good);
vecnum = vecnum(good);

for i = 1:length(vecfilenames)
    vecfilenames{i} = fullfile(basepath, vecdir, vecfilenames{i});
end

[vecnum,ord] = sort(vecnum);
vecfilenames = vecfilenames(ord);
imgfilenames = imgfilenames(ord);

fileskip = input('File skip (default = 1)? ','s');
if isempty(fileskip)
    fileskip = 1;
end

vecfilenames = vecfilenames(1:fileskip:end);
imgfilenames = imgfilenames(1:fileskip:end);
vecnum = vecnum(1:fileskip:end);

fprintf('Using files from %d to %d with frame skip %d\n', ...
    vecnum(1), vecnum(end), fileskip);

outdir = uigetdir(basepath, 'Choose output directory');

basename2 = input(sprintf('Use base output name ("%s")', basename), 's');
if ~isempty(basename2)
    basename = basename2;
end

timedWaitBar(0, 'Converting vector files...');
for i = 1:length(vecfilenames)
    [x,y,u,v,err,~,data] = loadInsight(vecfilenames{i});
    
    x = x/1000;     % convert to m
    y = y/1000;
    
    outvecfilename = sprintf('%s-vec-%06d.dat', basename, vecnum(i));
    dlmwrite(fullfile(outdir,outvecfilename), [x(:) y(:) u(:) v(:)], ',');
    timedWaitBar(i/length(vecfilenames));
end
timedWaitBar(1);

mperpix = data.MicrometersPerPixelX / 1e6;

x0 = data.OriginInImageX;
y0 = data.OriginInImageY;

fprintf('Delta T = %f sec\n', data.MicrosecondsPerDeltaT / 1e6);

trackLine(imgfilenames, [],[], 'savecallback',@saveoutline);


    function saveoutline(x,y, imname,imnum)
        
        outfilename = sprintf('%s-outline-%06d.dat', basename, vecnum(imnum));
        dlmwrite(fullfile(outdir,outfilename), [(x(:)+x0) (y(:)+y0)]*mperpix, ',');
        
    end

end
        