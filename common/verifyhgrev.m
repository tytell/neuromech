function diffs = verifyhgrev(hgold, varargin)

opt.datafilepath = '';
opt.comparehashes = false;
opt.quiet = false;

opt = parsevarargin(opt, varargin, 2);

analysisfiles = {};
csets = {};
revs = [];
maincsets = {};
mainrevs = [];
oldfiles = struct([]);
diffs = struct('type',{},'file',{},'old',{},'new',{});
for i = 1:length(hgold)
    if (isfield(hgold(i),'sub'))
        kv = opt2keyvalue(opt);
        diffs1 = verifyhgrev(hgold(i).sub, kv{:}, 'quiet');
        diffs = mergestruct(diffs, diffs1);
    end;
    
    oldfiles = mergestruct(oldfiles,hgold(i).datafiles);
    
    af1 = hgold(i).caller;
    cset1 = repmat({hgold(i).changeset},[1 length(af1)]);
    rev1 = repmat(hgold(i).rev, [1 length(af1)]);
    
    analysisfiles = [analysisfiles af1];
    csets = [csets cset1];
    revs = [revs rev1];
    if (isfield(hgold(i),'main'))
        maincsets = [maincsets hgold(i).main.changeset];
        mainrevs = [mainrevs hgold(i).rev];
    end;
end;

if (~isempty(opt.datafilepath))
    for i = 1:length(oldfiles),
        [~,fn,ext] = fileparts(oldfiles(i).name);
        oldfiles(i).name = fullfile(opt.datafilepath, [fn ext]);
    end;
end;

if (opt.comparehashes)    
    %compare data file hashes
    hgcurr = savehgrev([],'datafile',{oldfiles.name}, 'disablewarnings', 'hash',true);
    
    newfiles = hgcurr.datafiles;
    
    goodfiles = true(size(oldfiles));
    for i = 1:length(oldfiles)
        goodfiles(i) = strcmp(newfiles(i).sha1, oldfiles(i).sha1);
    end;
else
    %just compare data file dates and sizes
    hgcurr = savehgrev([],'datafile',{oldfiles.name}, 'disablewarnings', 'hash',false);
    
    newfiles = hgcurr.datafiles;
    
    goodfiles = (cat(2,newfiles.bytes) == cat(2,oldfiles.bytes)) & ...
        (cat(2,newfiles.datenum) == cat(2,oldfiles.datenum));
end;    

oldinfo = oldfiles(~goodfiles);
newinfo = newfiles(~goodfiles);

a = length(diffs)+1;
for i = 1:length(oldinfo),
    diffs(a).type = 'data';
    diffs(a).file = newinfo(i).name;
    diffs(a).old = oldinfo(i);
    diffs(a).new = newinfo(i);
    a = a+1;
end;
    
[ufiles,~,filecode] = unique(analysisfiles);
[ucsets,~,csetcode] = unique(csets);

[filecset,ind] = unique([filecode(:) csetcode(:)], 'rows');

a = length(diffs)+1;
for i = 1:size(filecset,1),
    diffs(a).type = 'analysis';
    diffs(a).file = ufiles{filecset(i,1)};
    diffs(a).old = {revs(ind(i)) ucsets{filecset(i,2)}};
    diffs(a).new = {hgcurr.rev hgcurr.changeset};
    a = a+1;
end;

isanalysis = cat(1,strcmp({diffs.type},'analysis'));

if (any(isanalysis))
    diffs1 = diffs(isanalysis);
    [ufiles,~,filecode] = unique({diffs1.file});
    
    csetsall = cell(length(diffs1),1);
    for i = 1:length(diffs1),
        csetsall{i} = diffs1(i).old{2};
    end;
    [ucsets,~,csetcode] = unique(csetsall);
    
    [~,ind] = unique([filecode(:) csetcode(:)], 'rows');
    diffs1 = diffs1(ind);

    if (any(~isanalysis))
        diffs = cat(1,diffs1(:), diffs(~isanalysis));
    else
        diffs = diffs1;
    end;
end;

if (~opt.quiet)
    if (isempty(diffs))
        fprintf('No differences!\n');
    else        
        if (strcmp(diffs(1).type, 'analysis'))
            fprintf('Analysis file differences.\n');
            
            i = 1;
            while ((i <= length(diffs)) && strcmp(diffs(i).type,'analysis')),
                [~,diff] = hg('diff -r ',[diffs(i).old{2} ':' hgcurr.changeset], ...
                    diffs(i).file, '--stat', 'echo',false);
                
                if (~isempty(diff)),
                    fprintf('  %s:\n    old %d:%s\n    new %d:%s\n', diffs(i).file, ...
                        diffs(i).old{1}, diffs(i).old{2}, diffs(i).new{1}, diffs(i).new{2});
                    fprintf('%s',diff);
                end;
                
                i = i+1;
            end;
        end;
    end;
end;
    