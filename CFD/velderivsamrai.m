function V = velderivsamrai(V)
% Calculates the velocity gradient derivatives in a samrai grid region V.
% Just a central difference.

N = numel(V);
for i = 1:N,
    [x1,y1] = meshgrid(V(i).px, V(i).py);

    V(i).dudx = NaN(size(x1));
    V(i).dudy = NaN(size(x1));
    V(i).dvdx = NaN(size(x1));
    V(i).dvdy = NaN(size(x1));

    V(i).dudx(:,2:end-1) = (V(i).U_0(:,3:end) - V(i).U_0(:,1:end-2)) ./ ...
        (x1(:,3:end) - x1(:,1:end-2));
    V(i).dudx(:,[1 end]) = (V(i).U_0(:,[2 end]) - V(i).U_0(:,[1 end-1])) ./ ...
        (x1(:,[2 end]) - x1(:,[1 end-1]));
    V(i).dudy(2:end-1,:) = (V(i).U_0(3:end,:) - V(i).U_0(1:end-2,:)) ./ ...
        (y1(3:end,:) - y1(1:end-2,:));
    V(i).dudy([1 end],:) = (V(i).U_0([2 end],:) - V(i).U_0([1 end-1],:)) ./ ...
        (y1([2 end],:) - y1([1 end-1],:));
    
    V(i).dvdx(:,2:end-1) = (V(i).U_1(:,3:end) - V(i).U_1(:,1:end-2)) ./ ...
        (x1(:,3:end) - x1(:,1:end-2));
    V(i).dvdx(:,[1 end]) = (V(i).U_1(:,[2 end]) - V(i).U_1(:,[1 end-1])) ./ ...
        (x1(:,[2 end]) - x1(:,[1 end-1]));
    V(i).dvdy(2:end-1,:) = (V(i).U_1(3:end,:) - V(i).U_1(1:end-2,:)) ./ ...
        (y1(3:end,:) - y1(1:end-2,:));
    V(i).dvdy([1 end],:) = (V(i).U_1([2 end],:) - V(i).U_1([1 end-1],:)) ./ ...
        (y1([2 end],:) - y1([1 end-1],:));
end;


