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

switch class(vals{1}),
    case 'double',
        addval = NaN;

    case 'char',
        addval = ' ';

    case 'logical',
        addval = false;
        
    case 'cell',
        addval = [];

    otherwise,
        error('Unsupported class');
end;


%run through and grow the original matrices
for a = 1:length(vals),
    v = vals{a};

    for d = 1:nd,
        if ((d ~= dim) && (size(v,d) < totalsz(d))),
            pmt = [d 1:d-1 d+1:nd];
            v = permute(v,pmt);
            if (iscell(v))
                sz1 = sz(pmt,a);
                v(end+1:totalsz(d),:) = cell(totalsz(d)-sz1(1),prod(sz1(2:end)));
            else
                v(end+1:totalsz(d),:) = addval;
            end;
            v = ipermute(v,pmt);
        end;
    end;
    
    vals{a} = v;
end;

%and do the concatenation
V = cat(dim, vals{:});

            
        

