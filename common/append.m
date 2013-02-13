function append(varargin)
% function append(file,variables...)
% Works like save, but appends to a file.
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell <tytell at jhu dot edu>

save(varargin{:},'-append');
