function [fn,ind] = sortfiles(fn)
% [fn,ind] = sortfiles(fn)
% Sorts files in a sensible way based on file numbering.
% demo1.tif, demo5.tif, demo10.tif comes out in that order, not as
% demo1, demo10, demo5, which is what sort would return.
% Note: This is probably rather slow for long lists of unnumbered
% files.

if (~isempty(fn)),
    [fn,ind] = sort(fn);
    for i = 1:length(fn),
        [rind,q,tok] = regexp(fn{i},'(.*[^0-9])([0-9]+)[^0-9]*');
        if (~isempty(rind)),
            tok = tok{1};
            base{i} = fn{i}(tok(1,1):tok(1,2));
            num(i) = str2num(fn{i}(tok(2,1):tok(2,2)));
        else
            base{i} = fn{i};
        end;
    end;

    [basenames,q,baseind] = unique(base);

    if (length(basenames) ~= length(fn)),
        for i = 1:length(basenames),
            k = find(baseind == i);
            if (length(k) > 1),
                [a,jig] = sort(num(k));
                fn(k) = fn(k(jig));
                ind(k) = ind(k(jig));
            end;
        end;
    end;
end;

       