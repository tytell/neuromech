function s = spaces(varargin)

if (length(varargin) == 1)
    sz = makerow(varargin{1});
elseif (all(cellfun(@numel,varargin) == 1))
    sz = cat(2,varargin{:});
else
    error('Invalid size vector');
end;
    
s = repmat(' ',sz);
