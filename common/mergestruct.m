function sm = mergeStruct(varargin)
% function sm = mergeStruct(a,b,c,...)
% Merges a,b,c, etc. into a structure array, with sm(1) being a, sm(2) being b,
% and so forth.  It assumes field names mostly overlap.

sm = varargin{1};

for i = 2:length(varargin),
	f1 = fieldnames(varargin{i});
	n = length(varargin{i});
	
	for j = 1:length(f1),
		sm = setfield(sm,{i:i+n-1},f1{j},getfield(varargin{i},f1{j}));
	end;
end;

			