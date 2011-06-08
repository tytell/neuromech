function [on,off,lab,len] = findruns(L,varargin)
%FINDRUNS  Finds runs of true values in a logical vector
% function [on,off,lab] = findruns(L,varargin)
%   Looks for repeated true values in a logical vector.  Returns
%   indices where the run started (on) and where it stopped (off).
%   Also produces a label vector, in which each run is labeled with an
%   increasing index.  Essentially a 1D version of BWLABEL.
%
% OPTIONS
%   'mingap' - Minimum length of a gap between runs to label them
%                    separately.  Default: 0
%
% SEE ALSO
%  BWLABEL

% Mercurial revision hash: $Revision$ $Date$
% See https://bitbucket.org/tytell/matlab/overview
% Copyright (c) 2011, Eric Tytell <tytell at jhu dot edu>

opt.mingap = 0;

opt = parsevarargin(opt,varargin,2);

if (size(L,1) == 1)
    istrans = true;
    L = L';
else
    istrans = false;
end;

on = find(~L(1:end-1) & L(2:end)) + 1;
off = find(L(1:end-1) & ~L(2:end));

if (L(1))
    on = [1; on];
end;
if (L(end))
    off = [off; length(L)];
end;

offlen = on(2:end) - off(1:end-1);

good = offlen > opt.mingap;
on = on([true; good]);
off = off([good; true]);

lab = zeros(size(L));
lab(on) = 1:length(on);
lab(off) = -(1:length(off));
lab = cumsum(lab);

if (L(end))
    lab(end) = lab(end-1);
end;

len = off - on + 1;

if (istrans)
    on = on';
    off = off';
    lab = lab';
    len = len';
end;

