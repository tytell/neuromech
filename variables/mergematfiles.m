function M = mergematfiles(infiles,vars,varargin)
%MERGEMATFILES   Merges variables in multiple .mat files
% M = mergematfiles(infiles,vars,...)
%
% Merges variables in multiple .mat files.  infiles specifies the files and vars specifies
% the variables to merge.
%
% Options:
%
%   'cat','catuneven','cell' - Specify one of these three to determine how the function
%   merges the variables.  It can use the function "cat", in which all of the variables
%   need to have the same size, or "catuneven" which produces a matrix with the size of
%   the largest individual variable, or it can put them all in a cell array.
%
%   'catvars','catunevenvars' - If you want to cat some vars and use catuneven on others,
%   you can specify the variable names after each of these options.
%
%   'missingvar' - What to do if a variable is missing.  'fail', 'skipvar', or 'skipfile'
%
%   'outfile' - Saves the output in the file specified.
%
%   'transposevectors' - How to handle vectors.  If this is true, it will transpose row
%   and column vectors so that they all align, and then merge them.  Otherwise, it will
%   merge them as is.
%
% Note: Supersedes FILECAT

% Mercurial revision hash: $Revision$ $Date$
% See https://bitbucket.org/tytell/matlab_variables/overview
% Copyright (c) 2011, Eric Tytell <tytell at jhu dot edu>

opt.mode = 'catuneven';
opt.catvars = {};
opt.catunevenvars = {};
opt.missingvar = 'fail';
opt.outfile = [];
opt.transposevectors = true;

opt = parsevarargin(opt,varargin,'multival',{'mode',{'cat','catuneven','cell'}});

%disable warnings for variables not found in a .mat file - we'll handle the problem
%ourselves
warnstate = warning('off','MATLAB:load:variableNotFound');

nvar = length(vars);
nfile = length(infiles);

mergevars = {};

emptystruct = cell(2,length(vars));
emptystruct(1,:) = vars;
A = cell(nvar,nfile);

%grab the variables from all of the files
for f = 1:nfile,
    F = load(infiles{f},vars{:});

    loadedvars = fieldnames(F);
    %check to make sure we got all of the variables
    if (length(loadedvars) ~= length(vars)),
        switch lower(opt.missingvar),
          case 'skipvar',
            missingvars = setdiff(vars,loadedvars);
            for i = 1:length(missingvars),
                F.(missingvars{i}) = [];
            end;
            
          case 'skipfile',
            F = struct(emptystruct{:});
            
          case 'fail',
            missingvars = setdiff(vars,loadedvars);
            missingvars = sprintf('%s ',missingvars{:});
            error('mergematfiles:badfile','File %s is missing vars: %s',...
                  infiles{f},missingvars);
        end;
    end;

    F = orderfields(F,vars);

    F = struct2cell(F);
    A(:,f) = F;
end;

M = struct(emptystruct{:});

%if they didn't specify specific variables for cat and catuneven, then use the mode to
%determine what to do
if (isempty(opt.catvars) && isempty(opt.catunevenvars)),
    switch lower(opt.mode),
      case 'cat',
        opt.catvars = vars;
        
      case 'catuneven',
        opt.catunevenvars = vars;
    end;
end;

%merge the variables
for i = 1:nvar,
    %check the dimensions
    nd = cellfun(@ndims,A(i,:));
    
    %and look for any variables that are empty (resulting from skipped files)
    good = ~cellfun(@isempty,A(i,:));
    firstgood = first(good);
    
    %transpose vectors if necessary
    if (opt.transposevectors && all(nd(good) == 2)),
        isrowvec = cellfun(@(x) ((numel(x) > 1) & (size(x,2) == numel(x))), A(i,:));
        iscolvec = cellfun(@(x) ((numel(x) > 1) & (size(x,1) == numel(x))), A(i,:));
        
        %if they're all vectors, but some are rows and some are columns, transpose them
        %so that they're all the same
        if (all(isrowvec(good) | iscolvec(good)) && ...
            any(isrowvec(good)) && any(iscolvec(good))),
            if (sum(isrowvec) > sum(iscolvec)),
                A(i,iscolvec) = cellfun(@transpose,A(i,iscolvec), ...
                                        'UniformOutput',false);
            else
                A(i,isrowvec) = cellfun(@transpose,A(i,isrowvec), ...
                                        'UniformOutput',false);
            end;
        end;
    end;
    
    %cat certain variables
    if (ismember(vars{i},opt.catvars)),
        iscat = false;
        nd1 = nd(firstgood);
        
        %check that they all have the same number of dimensions first
        if (all(nd(good) == nd1)),
            %and the same size
            sz = cellfun(@(x) (size(x)'), A(i,:), 'UniformOutput',false);
            sz = cat(2,sz{:});
            
            if (all(all(sz(:,good) == repmat(sz(:,firstgood),[1 sum(good)])))),
                A1 = cellfun(@(x) (shiftdim(x,-1)), A(i,:), ...
                             'UniformOutput',false);
                
                %merge only the good values along dimension 1 to begin with
                V0 = cat(1,A1{good});
                
                %build a matrix with space for good and bad values
                V1 = zeros([nfile sz(:,firstgood)'],class(A1{firstgood}));
                V1(good,:) = V0(:,:);
                
                iscat = true;
            end;
        end;
        
        if (~iscat),
            warning('mergematfiles:sizemismatch',...
                    'Size of var %s does not match in all files', vars{i});

            %if we couldn't cat the variables, we'll save them in a cell array
            V1 = A(i,:);
        end;
    elseif (ismember(vars{i},opt.catunevenvars)),
        %cat uneven
        A1 = cellfun(@(x) (shiftdim(x,-1)), A(i,:), ...
                     'UniformOutput',false);
        V0 = catuneven(1,A1{good});
        
        sz = size(V0);
        V1 = zeros([nfile sz(2:end)],class(A1{firstgood}));
        V1(good,:) = V0(:,:);
    else
        %save in a cell array
        V1 = A(i,:);
    end;

    %fill in missing values with NaNs or spaces, as appropriate
    switch class(V1),
      case 'double',
        V1(~good,:) = NaN;
      case 'char',
        V1(~good,:) = ' ';
    end;
    
    M.(vars{i}) = shiftdim(V1,1);
end;
        
%save the file names, too
M.files = infiles;
        
warning(warnstate);

if (~isempty(opt.outfile)),
    save(opt.outfile,'-struct','M');
end;

    