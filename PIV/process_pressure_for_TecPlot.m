function process_pressure_for_TecPlot(filenames, outlinefiles, outfilenames)

inputyn('', 'clearsaved',true);

if nargin == 0
    [fn,pathname] = uigetfile('*.dat', 'Choose pressure files', 'MultiSelect','on');
    filenames = cell(size(fn));
    for i = 1:length(fn)
        filenames{i} = fullfile(pathname,fn{i});
    end
    
    outlinefiles = regexprep(filenames, '-pressure-', '-outline-');
    outfilenames = regexprep(filenames, '-pressure-', '-tecplot-');
end

for i = 1:length(filenames)
    fprintf('%s\n', filenames{i});
    
    data = dlmread(filenames{i}, ',');

    d1 = find(diff(data(:,1)) ~= 0);
    d2 = length(d1)+1;

    if (d1(1)*d2 ~= size(data,1))
        shape = input(sprintf('What is the grid shape for the data? Appears to be about [%d, %d]: ',...
            d1(1), d2));
    else
        shape = [d1(1), d2];
    end

    if inputyn('Subtract median flow speed?', 'default',true)
        Umed = nanmedian(data(:,3));
        Vmed = nanmedian(data(:,4));
        
        fprintf('Umed = [%f, %f]\n', Umed,Vmed);
        
        data(:,3) = data(:,3) - Umed;
        data(:,4) = data(:,4) - Vmed;
    end
    
    if exist(outfilenames{i},'file') && ...
            ~inputyn('Overwrite file?', 'default',false)
        fprintf('Skipping!\n');
        continue;
    end
    
    outlinedata = dlmread(outlinefiles{i}, ',');
    
    x = data(:,1);
    y = data(:,2);
    ox = outlinedata(:,1);
    oy = outlinedata(:,2);
    x = reshape(x, shape(1), shape(2));
    y = reshape(y, shape(1), shape(2));
    
    [isfish,onedge] = inpolygon(x,y, ox,oy);
    isfish = isfish + onedge;
        
    fid = fopen(outfilenames{i},'w');

    fprintf(fid, 'Title = "%s"\n', filenames{i});
    fprintf(fid, 'Variables = X Y U V dpdx dpdy p absp bad isfish\n');
    fprintf(fid, 'ZONE I=%d, J=%d, F=POINT\n', shape(1), shape(2));

    bad = any(isnan(data),2);

    data(isnan(data)) = 0;
    data(:,9) = bad;

    data(:,10) = isfish(:);
    
    fprintf(fid, '%f, %f, %f, %f, %f, %f, %f, %f, %d, %d\n', data');

    fclose(fid);
end