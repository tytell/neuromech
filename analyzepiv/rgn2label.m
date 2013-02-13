function lab = rgn2label(x,y, rgns)
% function lab = rgn2label(x,y, rgns)
% Converts the Regions output of analyzePIV into a label matrix like
% bwlabeln would produce.

lab = zeros(size(x,1),size(x,2),size(rgns,2));
for j = 1:size(rgns,2),
    for i = 1:size(rgns,1),
        if (rgns(i,j).status > 0),
            in = inpolygon(x,y, rgns(i,j).x,rgns(i,j).y);

            % convert from a boolean to and index value, accounting
            % for the frame number
            in = find(in > 0);
            in = in + (j-1)*prod(size(x));
            lab(in) = i;
        end;
    end;
end;

    