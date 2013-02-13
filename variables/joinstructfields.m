function s = joinstructfields(s,app)
%JOINSTRUCTFIELDS   Joins the fields of one structure to another
%   S = joinstructfields(A,B)
%   Adds the fields of B (and their contents) to A.  Correctly handles
%   structure arrays, but A and B must be the same size.  If fields in A
%   are also present in B, then the values from B overwrite the values from
%   A.
%
% See also MAKESTRUCTARRAY

% Mercurial revision hash: $Revision$ $Date$
% See https://bitbucket.org/tytell/matlab_variables/overview
% Copyright (c) 2011, Eric Tytell <tytell at jhu dot edu>

if (isempty(s)),
    s = app;
elseif (isempty(app)),
    return;
else
    if ((ndims(s) ~= ndims(app)) || any(size(s) ~= size(app))),
        error('Structures to append must have the same size');
    end;

    fn = fieldnames(app);

    for i = 1:length(fn),
        [s.(fn{i})] = deal(app.(fn{i}));
    end;
end;
