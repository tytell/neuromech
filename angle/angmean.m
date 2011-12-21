function varargout = angmean(varargin)
% function [mean,R,stdev] = angmean(ang,dim)
%             or          = angmean(x,y,dim)
% second form is for unit vectors.  dim is optional, which makes the two
% forms ambiguous (in the rare instance that numel(ang) == 1).  To
% completely avoid ambiguity, you can use the form angmean(ang,[],dim).

if (nargin == 1),
  ang = varargin{1};
  angx = cos(ang);
  angy = sin(ang);
  dim = 1;
elseif (nargin == 2),
    if (numel(varargin{2}) == 1),
        [ang,dim] = deal(varargin{:});
        angx = cos(ang);
        angy = sin(ang);
    elseif ((ndims(varargin{1}) == ndims(varargin{2})) && ...
            all(size(varargin{1}) == size(varargin{2}))),
        [angx,angy] = deal(varargin{:});
        dim = 1;
    else
        error('x and y must have the same number of elements');
    end;
elseif (nargin == 3),
    [angx,angy,dim] = deal(varargin{:});
else
    error('Wrong number of parameters.');
end;

if (isempty(angy)),
    ang = angx;
    angx = cos(ang);
    angy = sin(ang);
end;

p = [dim 1:dim-1 dim+1:ndims(angx)];
angx = permute(angx,p);
angy = permute(angy,p);

angx = nanmean(angx);
angy = nanmean(angy);
a = atan2(angy,angx);
r = sqrt(angx.^2 + angy.^2);

s = sqrt(-2*log(r));

a = ipermute(a,p);
r = ipermute(r,p);
s = ipermute(s,p);

out = {a,r,s};
varargout = out(1:nargout);

