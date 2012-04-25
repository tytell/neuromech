function latextable(tab, varargin)
%LATEXTABLE  Produce a nice table in LaTeX format
%
%  latextable(table, options...)
%
% Produces a table in LaTeX format with sparklines and lots of options for spacing
% and formatting.  Takes a cell matrix or structure array containing the table
% items.  A cell matrix should be the table arranged as you would expect it
% with rows of the matrix becoming rows in the table.  In a structure, the
% structure fields are the columns and the elements in the structure array
% are the rows.  Since it uses LaTeX, any string item can have LaTeX
% formatting options (such as math operations).
%
% Sparklines can be defined as a vector of numbers, or in almost every
% format accepted by the plot command.  Note that the numbers are *not*
% scaled, and are expected to be in points.
%
% Many options are flexible in terms of size.  If an
% option has a single value, then it's applied to all the cells.  If it is
% in a row vector with the same number of columns as the table, then each
% separate item is applied to each column.  Same goes for a column vector
% with the same number of rows as the table.  And the option could also
% have the same number of elements as the whole table.
%
%OPTIONS
%   'rownames' or 'colnames' - Takes a cell array containing the names of the 
%      rows/columns.  Both can specify multicolumn/row labels.  See examples below.
%   'align' - 'c', 'l', or 'r', for center, left, or right.  Flexible size.
%   'numfmt','strfmt','multifmt' - Default formats for numbers, strings,
%      and cell elements in the table.  The 'format' option overrides these
%      values if it's specified.
%   'format' - printf style format specifier.  If empty, numbers, strings,
%      and cell items get the defaults specified in the options 'numfmt',
%      'strfmt', and 'multifmt'.  Sparklines are specified with 'spark'.
%      Formats can be complex for table entries with multiple values, such
%      as a format '%.1f (%.1f)' could be applied to a table cell with the
%      value {4.4, 0.2} to show values and errors, for instance.
%   'preamble','postamble' - LaTeX commands to apply before and after the
%     table.  Use at your own risk.
%   'tableenvironment' - LaTeX table environment.  Should be either
%     'longtable' (default) or 'table', unless you know LaTeX well.
%   'outfile' - Output file.  If not specified, output goes to the command
%     line.
%   'runlatex' - true or false to run the LaTeX command on the file.
%   'latex' - Path to latex executable ('/usr/texbin/pdflatex' by default)
%   'sparkwidth' - Default with of a sparkline in points.
%   'midrule' - Row numbers to put a line after.
%   'rowcolorstyle' - Empty or 'whitegray' for alternating white and gray
%     lines.  Or a cell array of valid LaTeX commands to define colors.
%     Will repeat as often as necessary.
%   'replacenan' - String with which to replace NaNs.
%   'header' - Text to put in the header on the page.
%   'orderfields' - For a structure table, names of the fields in the
%      desired order for the columns.  If empty, columns come in whatever
%      order the structure fields were defined.
%
%EXAMPLE
%  >> x = 1:10;  
%  >> y = x.^2;
%  >> y = y / max(y) * 15;      % scale y so that it's 15pts high
%  >> table = { 1 'A' {y}; 2 'B' {x,-y,'r:'}};
%  >> latextable(table);        % will show the table with default format
%                                 and no column or row names
%  >> latextable(table, 'colnames', {{2,'Stuff'},'','Spark'; ...
%     'Num','Str',''}, 'format',{'%03d','($%s^*$)','spark'})
%      % shows the table with a two row header, in which 'Stuff' spans the
%        first two columns.

opt.rownames = {};
opt.colnames = {};
opt.align = 'c';
opt.format = {};
opt.preamble = { ...
    '\documentclass[landscape,10pt]{article}' ...
    ' ' ...
    '\usepackage[margin=0.5in]{geometry}' ...
    '\usepackage{tikz,booktabs,longtable,multirow,colortbl}' ...
    '\renewcommand{\rmdefault}{ptm}' ...
    ' ' ...
    };
opt.postamble = { ...
    '\end{center}' ...
    '\end{document}' ...
    };
opt.tableenvironment = 'longtable';
opt.numfmt = '%g';
opt.strfmt = '%s';
opt.multifmt = 'spark';
opt.outfile = '';
opt.runlatex = true;
opt.showpdf = true;
opt.latex = '/usr/texbin/pdflatex';
opt.sparkwidth = 30;
opt.midrule = [];
opt.rowcolorstyle = 'whitegray';
opt.replacenan = '---';
opt.header = '';
opt.orderfields = {};

opt = parsevarargin(opt, varargin, 2, 'typecheck',false);

if (isnumeric(tab))
    tab = num2cell(tab);
elseif (isstruct(tab))
    n = numel(tab);
    if (~isempty(opt.orderfields))
        tab = getfieldsonly(tab, opt.orderfields);
    end
    tab = struct2cell(tab);
    tab = reshape(tab,[size(tab,1) n]);
    tab = tab';
    
    if (isstruct(opt.format))
        nfmt = numel(opt.format);
        if (~isempty(opt.orderfields))
            opt.format = getfieldsonly(opt.format, opt.orderfields);
        end
        opt.format = struct2cell(opt.format);
        opt.format = reshape(opt.format,[size(opt.format,1) nfmt]);
        opt.format = opt.format';
    end
end

ncol = size(tab,2);
nrow = size(tab,1);

if (ischar(opt.align) && (numel(opt.align) == 1))
    aligncell = opt.align(ones(1,ncol));
elseif (ischar(opt.align) && (length(opt.align) == ncol))
    aligncell = opt.align;
elseif (iscellstr(opt.align) && (length(opt.align) == ncol))
    aligncell = cat(2,opt.align{:});
end;

rownames = opt.rownames;
if (~isempty(rownames))
    if (numel(aligncell) == ncol)
        aligncell = ['l' aligncell];
    end;
    c0 = 1;
else
    c0 = 0;
end;

for i = 1:length(rownames),
    if (iscell(rownames{i}))
        if ((numel(rownames{i}) == 2) && ...
                isnumeric(rownames{i}{1}) && ischar(rownames{i}{2}))
            multirow = sprintf('\\multirow{%d}{*}{%s}', rownames{i}{:});
            
            rownames{i} = multirow;
        end;
    end;
end;

%add header if necessary
if (~isempty(opt.header))
    opt.preamble = [opt.preamble { ...
        '\usepackage{fancyhdr}' ...
        '\pagestyle{fancy}' ...
        ['\fancyhead[LE,LO]{' opt.header '}'] ...
        '\renewcommand{\footrulewidth}{0pt}' ...
        '\renewcommand{\headrulewidth}{0pt}' ...
        '' ...
        }];
end;
        
%sort out multicolumn column names
colnames = opt.colnames;
nhead = size(colnames,1);

if ((size(colnames,2) == ncol+1) && ~isempty(rownames))
    rowcolname = colnames(:,1);
    colnames = colnames(:,2:end);
else
    rowcolname = cell(nhead,1);
end;

cmidrule = cell(nhead-1,1);
isspan = false(nhead,ncol);
for i = 1:nhead,
    for j = 1:ncol,
        if (iscell(colnames{i,j}))
            if ((numel(colnames{i,j}) == 2) && ...
                    isnumeric(colnames{i,j}{1}) && ischar(colnames{i,j}{2}))
                nspan = colnames{i,j}{1};
                multicol = sprintf('\\multicolumn{%d}{c}{%s}', colnames{i,j}{:});
                
                colnames{i,j} = multicol;
                
                if (i < nhead)
                    gap = '';
                    if (j > 1)
                        gap = [gap 'l'];
                    end;
                    if (j+nspan-1 < ncol)
                        gap = [gap 'r'];
                    end;
                    isspan(i,j+(1:nspan-1)) = true;
                    
                    cmidrule{i} = [cmidrule{i} ...
                        sprintf('\\cmidrule(%s){%d-%d} ',gap,j+c0,j+c0+nspan-1)];
                end;
            else
                error('Unrecognized column name specification');
            end;
        end;
    end;
end;

fmt = cell(nrow,ncol);

numfmt = opt.numfmt;
strfmt = opt.strfmt;
multifmt = opt.multifmt;

if (ischar(opt.format))
    numfmt = opt.format;
elseif (iscell(opt.format) && (numel(opt.format) == ncol))
    fmt = repmat(opt.format(:)',[nrow 1]);
elseif (iscell(opt.format) && (numel(opt.format) == nrow))
    fmt = repmat(opt.format(:), [1 ncol]);
elseif (iscell(opt.format) && all(size(opt.format) == [nrow ncol]))
    fmt = opt.format;
elseif (iscell(opt.format) && (numel(opt.format) > 1))
    warning('latextable:fmtsize','Format is not empty but does not match table size.  Ignoring.');
end;

for i = 1:nrow,
    for j = 1:ncol,
        if (isnumeric(tab{i,j}))
            if (~isempty(opt.replacenan) && (numel(tab{i,j}) == 1) && isnan(tab{i,j}))
                tab{i,j} = opt.replacenan;
                fmt{i,j} = strfmt;
            elseif (isempty(fmt{i,j}))
                fmt{i,j} = numfmt;
            end;
        elseif (ischar(tab{i,j}))
            if (isempty(fmt{i,j}))
                fmt{i,j} = strfmt;
            end;
        elseif (iscell(tab{i,j}))
            if (isempty(fmt{i,j}))
                fmt{i,j} = multifmt;
            end;
        end;
    end;
end;

sparks = {};
a = 1;
for i = 1:nrow,
    for j = 1:ncol,
        if (isempty(tab{i,j}))
            fmt{i,j} = '';
            tab{i,j} = '';
        elseif (strcmpi(fmt{i,j},'spark'))
            if (isnumeric(tab{i,j}))
                plt = tab(i,j);
            else
                plt = tab{i,j};
            end;
            if ((numel(plt) == 1) || ...
                    (numel(plt) >= 2) && isnumeric(plt{1}) && ischar(plt{2}))
                y1 = plt{1};
                x1 = linspace(0,opt.sparkwidth,length(y1));
                p = 2;

                col = 'k';
                linestyle = '-';
                if (p <= numel(plt))
                    [islinespec,col1,mark,linestyle1] = matchlinespec(plt{p});
                    if (islinespec)
                        col = col1;
                        linestyle = linestyle1;
                        p = 3;
                    end;
                end;
                pltcmd = {x1,y1,col,linestyle};
            else
                p = 1;
                b = 1;
                pltcmd = cell(0,4);
                while ((p+1 <= numel(plt)) && isnumeric(plt{p}) && isnumeric(plt{p+1}))
                    [x1,y1] = deal(plt{p:p+1});
                    p = p+2;
                    
                    col = 'k';
                    linestyle = '-';
                    if (p <= numel(plt))
                        [islinespec,col1,mark,linestyle1] = matchlinespec(plt{p});
                        if (islinespec)
                            col = col1;
                            linestyle = linestyle1;
                            p = p+1;
                        end;
                    end;
                    pltcmd(b,:) = {x1,y1,col,linestyle};
                    b = b+1;
                end;
            end;
            
            nplt = size(pltcmd,1);
            
            pltopt.LineWidth = 0.75;
            pltopt = parsevarargin(pltopt, plt(p:end), p);
            
            if (~isempty(pltopt.LineWidth))
                drawspec0 = {sprintf('line width=%gpt',pltopt.LineWidth)};
            end;
            
            id = ['\spark' num2code(i,j)];
            sp1 = cell(4+nplt,1);
            sp1(1:2) = { ['\newcommand{' id '}{%'], ...
                '  \begin{tikzpicture} %' };
            
            for c = 1:nplt,
                drawspec = drawspec0;
                switch pltcmd{c,3}
                    case 'k',
                        % do nothing
                    case 'b',
                        drawspec = [drawspec {'blue'}];
                    case 'g',
                        drawspec = [drawspec {'green'}];
                    case 'r',
                        drawspec = [drawspec {'red'}];
                end;
                
                switch pltcmd{c,4}
                    case '-',
                        %do nothing
                    case '--',
                        drawspec = [drawspec {'dashed'}];
                    case ':',
                        drawspec = [drawspec {'dotted'}];
                end;
                
                x1 = pltcmd{c,1};
                y1 = pltcmd{c,2};
                if (~isempty(y1))
                    if (~isempty(drawspec))
                        str = '  \draw[';
                        for m = 1:length(drawspec)-1,
                            str = [str drawspec{m} ','];
                        end;
                        str = [str drawspec{end} '] '];
                    else
                        str = '  \draw ';
                    end;
                    
                    y1(~isfinite(y1)) = 0;
                    if (length(y1) > 1)
                        for m = 1:length(y1)-1,
                            str = [str sprintf('(%.2fpt,%.2fpt)--',x1(m),y1(m))];
                        end;
                    else
                        m = 0;
                    end;
                    str = [str sprintf('(%.2fpt,%.2fpt);%%',x1(m+1),y1(m+1))];
                    sp1{c+2} = str;
                else
                    sp1{c+2} = '';
                end;
            end;
            
            sp1{c+3} = '  \end{tikzpicture}%';
            sp1{c+4} = '}';
            
            sparks{a} = sp1;
            a = a+1;
            
            tab{i,j} = id;
            fmt{i,j} = '%s';
        end;
    end;
end;
            
if (isempty(opt.outfile))
    fid = 1;
else
    fid = fopen(opt.outfile,'w');
end;

fprintf(fid, '%s\n', opt.preamble{:});

if (~isempty(sparks))
    for i = 1:length(sparks)
        fprintf(fid, '%s\n', sparks{i}{:});
    end;
end;

begindoc = { '\begin{document}' ...
    '\begin{center}'};

fprintf(fid, '%s\n', begindoc{:});
fprintf(fid, '\\begin{%s}{%s}\n\\toprule\n\n', opt.tableenvironment, aligncell);

for i = 1:nhead,
    if (~isempty(rownames))
        fprintf(fid, '%s & ',rowcolname{i});
    end;
    for j = 1:ncol,
        if (~isspan(i,j))
            fprintf(fid, '%s', colnames{i,j});
            if ((j < ncol) && any(~isspan(i,j+1:end)))
                fprintf(fid, ' & ');
            else
                fprintf(fid, ' \\\\\n');
            end;
        end;
    end;
    if (i < nhead)
        fprintf(fid, '%s\n',cmidrule{i});
    end;
end;

fprintf(fid,'\\midrule \\endhead\n');


switch lower(opt.rowcolorstyle)
    case 'whitegray',
        rowcol = {'','\rowcolor[gray]{0.9}'};
    case 'none',
        rowcol = {''};
    otherwise,
        if (iscellstr(opt.rowcolorstyle))
            rowcol = opt.rowcolorstyle;
        else
            rowcol = {''};
        end;
end;
nrowcol = length(rowcol);

for i = 1:nrow,
    rc1 = rowcol{mod(i-1,nrowcol)+1};
    if (~isempty(rc1))
        fprintf(fid,'%s ',rc1);
    end;
    if (~isempty(rownames))
        fprintf(fid, '%s & ',rownames{i});
    end
    for j = 1:ncol,
        if (iscell(tab{i,j}))
            fprintf(fid, fmt{i,j}, tab{i,j}{:});
        else
            fprintf(fid, fmt{i,j}, tab{i,j});
        end;
        if (j < ncol)
            fprintf(fid, ' & ');
        else
            fprintf(fid, ' \\\\\n');
        end;
    end;
    if (ismember(i,opt.midrule))
        fprintf(fid,'\\midrule\n');
    end;
end;
fprintf(fid, '\\bottomrule\n');

fprintf(fid, '\\end{%s}\n', opt.tableenvironment);
fprintf(fid, '%s\n', opt.postamble{:});

if (fid ~= 1)
    fclose(fid);
    
    if (opt.runlatex)
        [pn,fn,ext] = fileparts(opt.outfile);
        if (isempty(pn))
            pn = '.';
        end;
        pdffile = fullfile(pn,[fn '.pdf']);

        if (exist(pdffile,'file'))
            delete(pdffile);
        end;
        
        done = false;
        rep = 1;
        while (~done && (rep < 5))
            [status,res] = system([opt.latex ' -halt-on-error -output-directory ' pn ...
                ' ' opt.outfile]);
            
            warnind = regexp(res,'Warning.*Rerun','once');
            if (isempty(warnind))
                done = true;
            elseif (status == 1)
                done = true;
            end;
            rep = rep+1;
        end;
        if (status == 1)
            error('LaTeX error');
        end;
        
        [pn,fn,ext] = fileparts(opt.outfile);
        pdffile = fullfile(pn,[fn '.pdf']);
        if (opt.showpdf && exist(pdffile,'file'))
            open(pdffile);
        end;
    end;
end;

        

function id = num2code(i,j)

ndig = max(1,ceil(log(i)/log(26)));
A = spaces(1,ndig);
for m = ndig:-1:1,
    A(m) = 'A' + mod(i-1,26);
    i = floor(i/26);
end;

ndig = max(1,ceil(log(j)/log(26)));
a = spaces(1,ndig);
for m = ndig:-1:1,
    a(m) = 'a' + mod(j-1,26);
    j = floor(j/26);
end;

id = [A a];
        

   



