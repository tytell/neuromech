function KV = opt2keyvalue(opt)
%OPT2KEYVALUE  Converts an opt structure to key value pairs
%
%  KV = opt2keyvalue(opt);
%
% See also PARSEVARARGIN

keynames = makerow(fieldnames(opt));
values = makerow(struct2cell(opt));

KV = [keynames; values];
KV = KV(:)';
