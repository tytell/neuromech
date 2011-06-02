function varargout = imshow6(varargin)
% h = imshow6(x,y,I)
%
% imshow with the version 6 parameter layout
%
% Mercurial revision hash: $Revision: 18f43cd9074e $ $Date: 2010/08/10 21:11:58 $
% Copyright (c) 2010, Eric Tytell

v = ver('images');
if (str2num(v.Version) >= 5),
    if ((nargin >= 3) && all(cellfun(@isnumeric,varargin(1:3)))),
        xd = varargin{1};
        yd = varargin{2};
        I = varargin{3};

        args = {I,'XData',xd,'YData',yd,varargin{4:end}};
        showaxes = true;
    else
        args = varargin;
        showaxes = false;
    end;

    chararg = find(cellfun(@ischar,args));
    notruesize = [strmatch('n',args(chararg),'exact') ...
        strmatch('notruesize',args(chararg),'exact')];

    if (numel(notruesize) == 1),
        notruesize = chararg(notruesize);
        args = args([1:notruesize-1 notruesize+1:end]);
        args = {args{:} 'InitialMagnification','fit'};
    end;

    h = imshow(args{:});
    if (showaxes),
        axis on;
    end;

else
    h = imshow(varargin{:});
end;
if (nargout == 1),
    varargout = {h};
end;

       