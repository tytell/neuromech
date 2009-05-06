function varargout = rowind(varargin)
% function ind = rowind(sz,r)
%   or     A = rowind(A,r)
%   or     [A,ind] = rowind(A,r)
%
% Returns an index that references a matrix of size sz at row r(i) and
% column i, for all columns i (or higher order dimensions).  Alternatively,
% you can pass a matrix A and it will return the value of the matrix at the 
% row values.   

if (nargin ~= 2)
    error('rowind:argnum','Wrong number of arguments');
end;

%parse the parameters
if ((numel(varargin{1}) == 2) && (numel(varargin{2}) == varargin{1}(2))),
    [sz,r] = deal(varargin{:});
elseif ((size(varargin{1},1) == 1) && ...
        (ndims(varargin{2}) == length(varargin{1})) && ...
        all(size(varargin{2}) == [1 varargin{1}(2:end)]))
    [sz,r] = deal(varargin{:});
elseif ((ndims(varargin{1}) == 2) && (numel(varargin{2}) == size(varargin{1},2))),
    [A,r] = deal(varargin{:});
    sz = size(A);
elseif ((ndims(varargin{1}) == ndims(varargin{2})) && ...
        (size(varargin{2},1) == 1) && ...
        (numel(varargin{2}) == numel(varargin{1})/size(varargin{1},1))),
    [A,r] = deal(varargin{:});
    sz = size(A);
else
    error('rowind:badargs','Unrecognized arguments');
end;

if (any((r < 0) | (r > sz(1)))),
    error('rowind:badrow','Row index out of range');
end;

%flatten higher order dimensions
sz0 = sz;
sz = [sz0(1) prod(sz0(2:end))];
r = reshape(r,[1 sz(2)]);

%make the index
ind = r + (0:sz(2)-1)*sz(1);

%return ind to the original shape
ind = reshape(ind,[1 sz0(2:end)]);

%and handle the output parameters
if (exist('A','var')),
    if (nargout == 2),
        varargout = {A(ind), ind};
    else
        varargout = {A(ind)};
    end;
else
    varargout = {ind};
end;
