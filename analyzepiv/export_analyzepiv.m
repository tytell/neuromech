function export_analyzepiv(filename, outname)

if nargin == 0
    [fn,pn] = uigetfile('*.mat','Select analyzepiv output file');
    filename = fullfile(pn,fn);
end

F = load(filename);
names = fieldnames(F);

names = names(~ismember(names,{'Regions','Units'}));
nfr = length(F.(names{1}));
len = zeros(1,length(names));
for i = 1:length(names)
    len1 = cellfun(@length,F.(names{i}));
    len1 = max(len1);
    len(i) = len1;
end

X = zeros(nfr,sum(len));
Xname = cell(1,sum(len));
a = 0;
for i = 1:length(names)
    good = ~cellfun(@isempty,F.(names{i}));
    
    x = cat(1,F.(names{i}){good});
    X(good,(1:len(i))+a) = x;
    
    if (len(i) == 1)
        Xname{a+1} = names{i};
    elseif (len(i) == 2)
        Xname{a+1} = [names{i} 'x'];
        Xname{a+2} = [names{i} 'y'];
    else
        for j = 1:len(i)
            Xname{a+j} = [names{i} num2str(j)];
        end
    end
    a = a+len(i);
end

tplt1 = repmat('%s, ',[1 sum(len)]);
tplt1 = [tplt1(1:end-2) '\n'];
tplt2 = repmat('%.4g, ',[1 sum(len)]);
tplt2 = [tplt2(1:end-2) '\n'];

[pn,fn,ext] = fileparts(filename);
if (nargin == 0)
    outname1 = fullfile(pn,[fn '.csv']);
    [fn,pn] = uiputfile('*.csv','Choose output file',outname1);
    outname = fullfile(pn,fn);
elseif (nargin == 1)
    outname = fullfile(pn,[fn '.csv']);
end

fid = fopen(outname,'w');
fprintf(fid,tplt1, Xname{:});
fprintf(fid,tplt2, X');
fclose(fid);



