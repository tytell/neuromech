function X = catfields(dim, S, fields)
%CATFIELDS  Concatenate fields of a structure
%
% Concatenates data from specified fields of a structure.

if (nargin == 3)
    isfld = ismember(fieldnames(S), fields);
else
    isfld = true(length(fieldnames(S)),1);
end

C = struct2cell(S);
sz = size(C);
C = C(isfld,:);

for i = 1:size(C,2)
    C{1,i} = cat(dim,C{:,i});
end

C = reshape(C(1,:),[1 sz(2:end)]);
X = cell2mat(C);