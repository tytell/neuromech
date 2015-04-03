function name = increment_file_name(basename,varargin)

opt.maxscan = 100;
opt = parsevarargin(opt,varargin, 2);

[pn,fn,ext] = fileparts(basename);
[numstr,a,b] = regexp(fn,'\d+','match','start','end');

if (~isempty(numstr))
    numstr = numstr{end};
    a = a(end);
    b = b(end);
    len = b-a+1;
    
    num = str2double(numstr);
else
    num = 1;
    len = 3;
    a = length(fn)+1;
    b = length(fn);
end

done = false;
i = 0;
while ~done
    fn1 = sprintf('%s%0*d%s%s',fn(1:a-1),len,num+i,fn(b+1:end),ext);
    name = fullfile(pn,fn1);
    
    if ~exist(name, 'file')
        done = true;
    else
        i = i+1;
        if (i > opt.maxscan)
            error('Too many files');
        end
    end
end
