function tf = inputyn(str,default)
% function tf = inputyn(str,default)
% Returns true if the user types anything with a 'y' or 'Y' in it
% Returns default if the user doesn't type anything.
% Default is optional and defaults to true.
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell

if (nargin == 1),
	default = 1;
end;

yn = input(str,'s');
if (isempty(yn)),
	tf = default;
else
	tf = ~isempty(strmatch(yn,{'y','Y'}));
end;

