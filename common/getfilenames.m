function filenames = getfilenames(d,varargin)
% function filenames = getfilenames(d,options...)
% Returns filenames on path d (potentially with wildcards).
% 
% Options:
%   'recursive' - Recurses through directories
%   'exclude' - Cell array of exact file patterns (not wildcards) that
%     should be excluded from results.  Default is {'.','..','.DS_Store'}
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell <tytell at jhu dot edu>

opt.recursive = false;
opt.exclude = {};
opt.baseexclude = {'^\.\.?$','^\.DS_Store$','^\.LINK'};
opt.include = {};

opt = parsevarargin(opt,varargin,2);

if (exist(d,'dir'))
    pn = d;
    fn = '';
    ext = '';
else
    [pn,fn,ext] = fileparts(d);
end;

exclude = [opt.exclude opt.baseexclude];

files = dir(d);
if (~isempty(files)),
    a = 1;
    filenames = cell(length(files),1);
	for i = 1:length(files),
        fn1 = fullfile(pn,files(i).name);
        good = true;
        for j = 1:length(exclude),
            if (~isempty(regexp(files(i).name, exclude{j}, 'once')))
                good = false;
                break;
            end;
        end;
        if (good && ~isempty(opt.include))
            if (opt.recursive && exist(fn1,'dir'))
                good = true;
            else
                good = false;
                for j = 1:length(opt.include)
                    if (~isempty(regexp(files(i).name, opt.include{j},'once')))
                        good = true;
                        break;
                    end
                end
            end
        end
        if (good)
            if (opt.recursive && exist(fn1,'dir')),
                rec = getfilenames(fullfile(fn1,[fn ext]),'recursive',opt.recursive, ...
                    'exclude',opt.exclude, 'baseexclude',opt.baseexclude);
                if (~isempty(rec))
                    filenames(a:a+length(rec)-1,1) = rec;
                end;
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
