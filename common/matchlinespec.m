function [islinespec,col,mark,linestyle] = matchlinespec(str)

[col,split] = regexp(str, '[rgbcmykw]', 'once','match','split');
left = strcat(split{:});

[linestyle,split] = regexp(left, '(--|-.|-|:)', 'once','match','split');
left = strcat(split{:});

[mark,split] = regexp(left, '(square|diamond|pentagram|hexagram|[+o*.xsd^v><ph])', ...
    'once','match','split');
left = strcat(split{:});

islinespec = isempty(left);
