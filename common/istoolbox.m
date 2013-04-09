function t = istoolbox(name,version)

if nargin == 1
    version = 0;
end

v = ver;
ind = find(~cellfun(@isempty,regexpi({v.Name}, name)));
if (length(ind) == 1)
    vinst = str2double(regexp(v(ind).Version, '\.','split'));
    
    if (isnumeric(version))
        version = num2str(version);
    end
    vreq = str2double(regexp(version,'\.','split'));
    
    if (length(vreq) > length(vinst))
        vinst(length(vreq)) = 0;        % pad with zeros
    elseif (length(vinst) > length(vreq))
        vreq(length(vinst)) = 0;
    end
    
    t = all(vinst >= vreq);
else
    t = false;
end
