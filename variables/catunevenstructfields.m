function S = catunevenstructfields(dim,varargin)
%CATSTRUCTFIELDS   Concatenates matching fields of multiple structures

if (nargin == 0)
    S = struct([]);
else
    S = varargin{1};
    for i = 2:length(varargin)
        S1 = varargin{i};
        fn = fieldnames(S1);
    
        for i = 1:length(fn)
            if (~isfield(S, fn{i}))
                error('Fields of structures do not match');
            end
            [S.(fn{i})] = catuneven(dim,S.(fn{i}),S1.(fn{i}));
        end
    end
end
