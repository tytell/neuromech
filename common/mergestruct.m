function sm = mergestruct(varargin)
% function sm = mergeStruct(a,b,c,...)
% Merges a,b,c, etc. into a structure array, with sm(1) being a, sm(2) being b,
% and so forth.  It assumes field names mostly overlap.
%
% Option 'flat' just merges all of the fields into a single structure.
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell <tytell at jhu dot edu>

isflat = false;
if ((nargin > 1) && ischar(varargin{end})),
    switch lower(varargin{end}),
        case 'flat',
            isflat = true;
        otherwise,
            error('Unrecognized option %s.', varargin{end})
    end;
    varargin = varargin(1:end-1);
end;

sm = varargin{1};

if (isflat),
    k = 1;
else
    k = length(sm)+1;
end;

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
    if (~isflat),
        k = k+n;
    end;
end;

			