function sm = makestructarray(varargin)
%MAKESTRUCTARRAY   Merges similar (but not identical) structures into an array
%  sm = makestructarray(a,b,c,...)
%  Merges a,b,c, etc. into a structure array, with sm(1) being a, sm(2) being b,
%  and so forth.  It assumes field names mostly overlap.
%
%  Option 'flat' is no longer allowed.  Use joinstructfields instead.
%
% See also JOINSTRUCTFIELDS.

% Mercurial revision hash: $Revision$ $Date$
% See https://bitbucket.org/tytell/matlab_variables/overview
% Copyright (c) 2011, Eric Tytell <tytell at jhu dot edu>

sm = varargin{1};

k = length(sm)+1;

for i = 2:length(varargin),
	f1 = fieldnames(varargin{i});
	n = length(varargin{i});

    if (isempty(fieldnames(sm)) && isempty(f1)),
        sm(k:k+n-1) = struct;
    elseif (isempty(f1)),
        f1 = fieldnames(sm);
        sm(k:k+n-1).(f1{1}) = deal([]);
    else
        for j = 1:length(f1),
            [sm(k:k+n-1).(f1{j})] = deal(varargin{i}.(f1{j}));
        end;
    end;
    k = k+n;
end;

			