function C = camcopy(varargin)
% function camcopy(from,to, options to copy...)
% or
%  C = camcopy('from',figure(2));
%  camcopy('to',figure(3),C);
%
% Copies camera options from one figure to another.  By default, copies all
% camera-related options ('CameraPosition','CameraPositionMode', 'CameraTarget',
% 'CameraTargetMode','CameraUpVector','CameraUpVectorMode','CameraViewAngle',
% 'CameraViewAngleMode', and 'Projection') so that two figures can have
% identical 3D views.  Can also select which options to copy:
%
% C = camcopy('from',fig1);
% camcopy('to',fig2,C,'CameraPosition','CameraTarget');
%  -> only copies 'CameraPosition' and 'CameraTarget'
%
% C = camcopy('from',fig1);
% camcopy('to',fig2,C,'except','Projection');
%  -> copies all options except for 'Projection'
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell <tytell at jhu dot edu>



allopts = {'CameraPosition','CameraPositionMode', ...
           'CameraTarget','CameraTargetMode',...
           'CameraUpVector','CameraUpVectorMode',...
           'CameraViewAngle','CameraViewAngleMode',...
           'Projection'};

from = [];
to = [];
if (ischar(varargin{1})),
    done = false;
    i = 1;
    while ((i <= nargin) && ~done && ischar(varargin{i})),
        switch lower(varargin{i}),
          case 'from',
            from = varargin{i+1};
            i = i+2;
            
          case 'to',
            to = varargin{i+1};
            i = i+2;
            
          otherwise,
            done = true;
        end;
    end;
else
    if (nargin < 2),
        error('Too few arguments');
    end;
    from = varargin{1};
    to = varargin{2};
    i = 3;
end;
    
opts = varargin(i:end);
if (isempty(opts)),
    opts = allopts;
elseif (isstruct(opts{1})),
    vals = struct2cell(opts{1});
    opts1 = fieldnames(opts{1});
    opts = [opts1 opts(2:end)];
end;

ind = strmatch('except',opts,'exact');
if (~isempty(ind)),
    opts = varargin(1:ind-1);
    exceptopts = varargin(ind+1:end);
else
    exceptopts = {}; 
end;

opts = opts(~ismember(opts,exceptopts));

if (~isempty(from) && strcmp(get(from,'Type'),'figure')),
    ax = findobj(from,'Type','axes');
    if (length(ax) > 1)
        error('Cannot figure out which axes you mean');
    end;
    from = ax;
end;
if (~isempty(to) && strcmp(get(to,'Type'),'figure')),
    to = findobj(to,'Type','axes');
end;

for i = 1:length(opts),
    if (~isempty(from)),
        val = get(from,opts{i});
    else
        val = vals{i};
    end;
    if (~isempty(to)),
        set(to,opts{i},val);
    end;
    
    C.(opts{i}) = val;
end;

    