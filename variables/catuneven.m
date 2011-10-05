function V = catuneven(dim, varargin)
%CATUNEVEN  Concatenates matrices when some dimensions may not match
%  V = catuneven(dim, V1, V2, ...)
%  Concatenates matrices along the dimension dim, when the other dimensions
%  may not match in size.  Grows the final matrix so that it's filled with
%  NaNs where the original matrices were too small.
%
% See also CAT.

% Mercurial revision hash: $Revision$ $Date$
% See https://bitbucket.org/tytell/matlab_variables/overview
% Copyright (c) 2011, Eric Tytell <tytell at jhu dot edu>

%nothing to cat?  Then return an empty matrix
if (nargin == 1),
    V = [];
    return;
end;

vals = varargin;
%get rid of empty matrices
% good = ~cellfun(@isempty,vals);
% vals = vals(good);

%if we got rid of all of the matrices, or left only one, deal with them
%specially
if (isempty(vals)),
    V = [];
    return;
elseif (length(vals) == 1),
    V = vals{1};
    return;
end;

%get the maximum number of dimensions
nd = cellfun(@ndims, vals);
nd = nd(1);
if (nd < dim)
    nd = dim;
end;

%save the sizes of everything
sz = ones(nd, length(vals));
for d = 1:nd,
    sz(d,:) = cellfun('size', vals, d);
end;

%and figure out the final output size
totalsz = max(sz,[],2);
totalsz(dim) = sum(sz(dim,:));
totalsz = totalsz';

switch class(vals{1}),
    case 'double',
        addval = NaN;

    case 'char',
        addval = ' ';

    case 'logical',
        addval = false;
        
    case 'cell',
        addval = [];
        
    case {'int8','uint8','int16','uint16','int32','uint32', ...
            'int64','uint64'}
        addval = zeros(1,1,class(vals{1}));

    otherwise,
        error('Unsupported class');
end;

totalsz1 = totalsz;
totalsz1(dim) = 1;
switch class(vals{1}),
    case 'double',
        vs = NaN(totalsz1);
    case 'char',
        vs = repmat(' ',totalsz1);
    case 'logical',
        vs = false(totalsz1);
    case 'cell',
        vs = cell(totalsz1);
    case {'int8','uint8','int16','uint16','int32','uint32', ...
            'int64','uint64'}
        vs = zeros(totalsz1,class(vals{1}));
        
    otherwise,
        error('Unsupported class');
end;

%run through and grow the original matrices
for a = 1:length(vals),
    v = vals{a};
    sz1 = size(v);
    
    if (~isempty(v))
        repdim = ones(1,nd);
        repdim(dim) = size(v,dim);
        
        vs1 = repmat(vs, repdim);
        
        %good sets the area within vs1 in which v fits
        good = true(totalsz1);
        for d = 1:nd,
            %bring the d dimension first
            pmt = [d 1:d-1 d+1:nd];
            good = permute(good, pmt);
            good1 = false(size(good));
            if (d ~= dim)
                good1(1:sz1(d),:) = true;
            else
                good1(1,:) = true;
            end;
            %make sure we are true along other dimensions, too
            good = good & good1;
            good = ipermute(good, pmt);
        end;
        good = repmat(good,repdim);
        
        vs1(good) = v;
        
        vals{a} = vs1;
    end;
end;

%and do the concatenation
V = cat(dim, vals{:});

            
        

