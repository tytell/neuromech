function sm = makestructarray(varargin)
%MAKESTRUCTARRAY   Merges similar (but not identical) structures into an array
%  sm = makestructarray(a,b,c,...)
%  Merges a,b,c, etc. into a structure array, with sm(1) being a, sm(2) being b,
%  and so forth.  It assumes field names mostly overlap.
%
%  Option 'skipempty' defines how empty structures are handled.  Default
%  (skipempty = false) fills empty structure fields so that the length of
%  the final array is equal to the number of structures passed in.
%  
%    e.g. >> S1 = struct('test',4,'thing',12);
%         >> S = makestructarray(S1,[],S1)
%         S = 
%         1x3 struct array with fields:
%             test
%             thing
%         >> S(2)
%         ans = 
%              test: []
%              thing: []
%
% but
%         >> S = makestructarray(S1,[],S1,'skipempty')
%         S = 
%         1x2 struct array with fields:
%             test
%             thing
% (skipempty=true is similar to the behavior of cat with empty elements)
%
%  Option 'flat' is no longer allowed.  Use joinstructfields instead.
%
% See also JOINSTRUCTFIELDS.

% Mercurial revision hash: $Revision$ $Date$
% See https://bitbucket.org/tytell/matlab_variables/overview
% Copyright (c) 2011, Eric Tytell <tytell at jhu dot edu>

opt.skipempty = false;

[opt,S] = parsevarargin(opt,varargin,'leaveunknown');

if (length(S) >= 1)
    sm = S{1};
    
    k = numel(sm)+1;
    
    for i = 2:length(S),
        if (isempty(S{i}))
            if (opt.skipempty)
                n = 0;
                f1 = {};
            else
                n = 1;
                f1 = {};
            end;
        else
            f1 = fieldnames(S{i});
            n = numel(S{i});
        end;
        
        if (n > 0)
            if (isempty(fieldnames(sm)) && isempty(f1)),
                sm(k:k+n-1) = struct;
            elseif (isempty(f1)),
                f1 = fieldnames(sm);
                sm(k:k+n-1).(f1{1}) = deal([]);
            else
                for m = 1:n,
                    sm(k+m-1).(f1{1}) = S{i}(m).(f1{1});
                end
                for j = 2:length(f1),
                    [sm(k:k+n-1).(f1{j})] = deal(S{i}.(f1{j}));
                end;
            end;
            k = k+n;
        end;
    end;
else
    sm = struct({});
end;

			