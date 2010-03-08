function filenames = getfilenames(d,varargin)

opt.recursive = false;
opt.exclude = {'.','..','.DS_Store'};

opt = parsevarargin(opt,varargin,2);

if (exist(d,'dir'))
    pn = d;
    fn = '';
    ext = '';
else
    [pn,fn,ext] = fileparts(d);
end;

files = dir(d);
if (~isempty(files)),
    a = 1;
    filenames = cell(length(files),1);
	for i = 1:length(files),
        if (~ismember(files(i).name, opt.exclude)),
            fn1 = fullfile(pn,files(i).name);
            if (opt.recursive && exist(fn1,'dir')),
                rec = getfilenames(fullfile(fn1,[fn ext]));
                filenames(a:a+length(rec)-1,1) = rec;
                a = a+length(rec);
            else
                filenames{a,1} = fullfile(pn,files(i).name);
                a = a+1;
            end;
        end;
	end;
    filenames = filenames(1:a-1);
else
	filenames = {};
end;

filenames = sortfiles(filenames);
