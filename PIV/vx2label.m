function varargout = vx2Label(vxr,x,y,nfr)
% function L = vx2Label(vxr,x,y,nfr)
% or       [L,names] = vx2Label(vxr,x,y,nfr)

if (isfield(vxr,'name')),
    names = {vxr.name};
    [names,q,nameind] = unique(names);
    
    %deal with empty elements.  They will always be nameind = 1, because
    %of how unique sorts them.  We want them all to be at the end of the
    %named vortices
    nameind = nameind-1;
    emptyind = length(names);

    k = find(nameind == 0);
    nameind(k) = emptyind;
    names = names(2:end);
else
    nameind = 1:length(vxr);
end;

nameind = uint8(nameind);
L = uint8(zeros(size(x,1),size(x,2),nfr));
szxy = prod(size(x));

phi = linspace(0,2*pi,17)';
x1 = zeros(length(phi),1);
y1 = zeros(length(phi),1);

timedWaitBar(0, 'Generating label matrix...');
for i = 1:length(vxr),
    nfr = length(vxr(i).frames);

    %get the angle mod pi
    ang1 = mod(vxr(i).angle,pi);

    %this is the value that phi will have when we're actually at the
    %angle ang1.  It's different because a ~= b
    a = vxr(i).majordiam/2;
    b = vxr(i).minordiam/2;
    ang2 = atan2(a.*sin(ang1), b.*cos(ang1));

    for j = 1:nfr,
        x1 = a(j)*cos(phi-ang2(j)).*cos(ang1(j)) - ...
             b(j)*sin(phi-ang2(j)).*sin(ang1(j)) + ...
             vxr(i).ctrx(j);
        y1 = a(j)*cos(phi-ang2(j)).*sin(ang1(j)) + ...
             b(j)*sin(phi-ang2(j)).*cos(ang1(j)) + ...
             vxr(i).ctry(j);

        %find elements inside the outline (x1,y1)
        in = inpolygon(x,y, x1,y1);
        in = find(in > 0);
        
        %shift them up to be on the right page of L
        in = in + (vxr(i).frames(j) - 1)*szxy;
        L(in) = nameind(i);
    end;

    timedWaitBar(i/length(vxr), 'Generating label matrix...');
end;
timedWaitBar(1);

varargout{1} = L;
if (exist('names') & (nargout == 2)),
    varargout{2} = names;
end;
