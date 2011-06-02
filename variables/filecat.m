function filecat(dim,files, varargin)
%FILECAT   Concatenate variables with the same name in a set of files
%   filecat(dim,files)
%   Looks for variables with the same name in each file in files and
%   concatenates them along the dimension dim using the same procedure as
%   CATUNEVEN.  Saves the data in the main workspace.
%
% Note: Better to use MERGEMATFILES.
%
% See also CATUNEVEN and MERGEMATFILES.

% Mercurial revision hash: $Revision$ $Date$
% See https://bitbucket.org/tytell/matlab_variables/overview
% Copyright (c) 2011, Eric Tytell <tytell at jhu dot edu>


%load the first file to get the variable names and sizes
F = load(files{1});
C = struct2cell(F);
varnames = fieldnames(F);
nvar = length(C);

%collect the variable sizes
ndmax = max(cellfun(@ndims,C));
ndmax = max(ndmax,dim);
for i = 1:ndmax,
    sz(:,i) = cellfun('size',C,i)';
end;
%index for dimensions that aren't the one we're interested in
otherdims = [1:dim-1 dim+1:ndmax];

for fn = 2:length(files),
    %now load the next files
    F1 = load(files{fn});
    C1 = struct2cell(F1);
    varnames1 = fieldnames(F1);
    
    for var = 1:length(varnames1),
        %look for variables with the same name
        match = strmatch(varnames1{var},varnames,'exact');
        if (~isempty(match)),
            sz1 = size(C1{var});
            sz1(end+1:ndmax) = 1;
            
            %check to see whether the variables are numeric and if any of
            %the sizes match
            if (isnumeric(C1{var}) && isnumeric(C{match}) && ...
                    any((sz1 > 1) & (sz1 == sz(match,:)))),
                
                %expand some dimensions if necessary
                for j = 1:ndmax-1,
                    j = otherdims(j);
                    
                    short = size(C1{var},j) - sz(match,j);
                    if (short > 0),
                        sznew = sz(match,:);
                        sznew(j) = abs(short);
                        C{match} = cat(j,C{match},repmat(NaN,sznew));
                        sz(match,j) = size(C{match},j);
                    elseif (short < 0),
                        sznew = size(C1{var});
                        sznew(j) = abs(short);
                        C1{var} = cat(j,C1{var},repmat(NaN,sznew));
                    end;
                end;
                %cat the variables
                C{match} = cat(dim, C{match}, C1{var});
            elseif (isa(C{match}, class(C1{var})) && ...
                    all(sz1(otherdims) == sz(match,otherdims))),
                %for other classes, just cat them normally
                C{match} = cat(dim, C{match}, C1{var});
            else
                %otherwise, make a cell array of the variables
                if (fn == 2),
                    C{match} = {C{match}};
                end;
                C{match}{fn} = C1{var};
            end;               
        else
            %if we didn't find a match for the variable name, add it in
            nvar = nvar+1;
            C{nvar} = C1{var};
            varnames{nvar} = varnames1{var};
            sz1 = size(C1{var});
            sz1(end+1:ndmax) = 1;
            sz(nvar,:) = sz1;
        end;
    end;
end;

%save all the variables in the main workspace
for var = 1:length(varnames),
    assignin('base',varnames{var},C{var});
end;

