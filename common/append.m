function append(varargin)
% function append(file,variables...)
% Works like save, but appends to a file.
%
% Mercurial revision hash: $Revision: 18f43cd9074e $ $Date: 2010/08/10 21:11:58 $
% Copyright (c) 2010, Eric Tytell <tytell at jhu dot edu>

save(varargin{:},'-append');
