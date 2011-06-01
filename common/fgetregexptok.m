function varargout = fgetregexptok(fid, expr, ln)

if (nargin == 2)
    ln = '';
end;

if (isempty(ln))
    ln = fgetl(fid);
end;

tok = regexp(ln, expr, 'tokens','once');
if (nargout == length(tok)+1)
    varargout = [tok ''];
elseif (nargout == length(tok))
    varargout = tok;
else
    varargout = [cell(1,nargout-1) ln];
end;


