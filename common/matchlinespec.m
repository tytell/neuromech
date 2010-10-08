function [islinespec,col,mark,linestyle] = matchlinespec(str)
% function [islinespec,col,mark,linestyle] = matchlinespec(str)
% Tests whether a string matches a Matlab style linespec.  Returns true in
% islinespec if the linespec is valid.  Also returns the separated color,
% marker, and linestyle.
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell <tytell at jhu dot edu>

[col,split] = regexp(str, '[rgbcmykw]', 'once','match','split');
left = strcat(split{:});

[linestyle,split] = regexp(left, '(--|-.|-|:)', 'once','match','split');
left = strcat(split{:});

[mark,split] = regexp(left, '(square|diamond|pentagram|hexagram|[+o*.xsd^v><ph])', ...
    'once','match','split');
left = strcat(split{:});

islinespec = isempty(left);
