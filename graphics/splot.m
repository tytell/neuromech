function varargout = splot(varargin)
% function h = splot(...)
%
% Just like plot, but runs the function squeeze on each of its arguments
% before it calls plot.
%
% Useful for situations when you have a high dimensional data set.  For
% example, if size(M) = [20 10 3 4] and size(N) = [20 1 3 4], then
% 
% splot(M(:,3,:,2),N(:,1,:,2))
%
% works in the way it seems it should, whereas plot would complain.

% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2011, Eric Tytell <tytell at jhu dot edu>


v = cell(size(varargin));
for i = 1:length(varargin),
	v{i} = squeeze(varargin{i});
end;

h = plot(v{:});
if (nargout == 1)
	varargout{1} = h;
end;

