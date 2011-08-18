function [S,good] = getfieldsonly(S, varargin)
%GETFIELDSONLY    Gets a structure with only certain fields
%   S = getfieldsonly(S, fields, ...)
%   Gets a structure that contains only the fields specified in fields

% Mercurial revision hash: $Revision$ $Date$
% See https://bitbucket.org/tytell/matlab_variables/overview
% Copyright (c) 2011, Eric Tytell <tytell at jhu dot edu>

if ((nargin == 2) && iscellstr(varargin{1})),
    fields = varargin{1};
else
    fields = varargin;
end;

fn = fieldnames(S);

good = ismember(fields,fn);
if (any(~good) && (nargout ~= 2))
    warning('getfieldsonly:missingfields','Not all requested fields are present in the structure');
end;
fields = fields(good);

fn = setdiff(fn, fields);

S = rmfield(S, fn);
S = orderfields(S,fields);

good = all(good);



