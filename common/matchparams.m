function tp = matchparams(p, req, opt)

nd = 0;
for i = 1:length(p),
    nd = max(nd,ndims(p{i}));
    sz(1:ndims(p{i}),i) = size(p{i})';
end;

for i = 1:length(req),
    r1 = req{i};

    [match(i),len(i),szlist,r1] = matchparamtype(p,sz,req{i}, []);
    req{i} = r1;
end;

if (sum(match) == 0),
    tp = 0;
elseif (sum(match) == 1),
    tp = find(match == 1);
else
    [q,tp] = max(len.*match);
end;

r1 = req{tp};
for i = 1:length(r1),
    pmatch = r1{i}{end};
    if (isnumeric(pmatch) & (length(pmatch) == 1)),
        v = p{pmatch};
    else
        for j = 1:length(pmatch),
            if (length(pmatch{j}) == 1),
                v{j} = p{pmatch{j}};
            else
                v{j} = p(pmatch{j});
            end;
        end;
    end;
    assignin('caller',r1{i}{1},v);
end;

% -------------------------------------
function [match,len,szlist,r1] = matchparamtype(p,sz, r1, szlist)

match = 1;
len = 0;
for j = 1:length(r1),
    pspec = r1{j};

    if (iscell(pspec{2})),
        k = 1;
        m0 = 0;
        pspec2 = pspec;
        while (~m0 & (k <= length(pspec{2}))),
            pspec2{2} = pspec{2}{k};
            [m0,l,szlist] = matchparamtype(p(j),sz(:,j), {pspec2}, ...
                                           szlist);
            k = k+1;
        end;
        if (~m0),
            match = 0;
            break;
        end;
        pspec{end+1} = j;
    else
        done = 0;
        switch pspec{2},
         case {'*', '+'},
          if (iscell(pspec{3}) & (length(pspec) == 3)),
              r2 = pspec{3};
          else
              r2 = {pspec([1 3:end])};
          end;

          [m0,l,szlist] = matchparamtype(p(j:end),sz(:,j:end), r2, ...
                                         szlist);
          if (~m0 & (pspec{2} == '+')),
              match = 0;
              break;
          end;
          pspec{end+1} = {j:j+l-1};

          k = j+l;
          while (m0 & (k <= length(p))),
              [m0,l,szlist] = matchparamtype(p(k:end),sz(:,k:end), ...
                                             r2, szlist);
              pspec{end} = {pspec{end} k:k+l-1};
              k = k+l;
          end;
          len = k-1;
          done = 1;

         case 'cellstr',
          if (~iscellstr(p{j})),
              match = 0;
              break;
          end;
          matchsz = pspec{3};

         case 'cell',
          if (~iscell(p{j})),
              match = 0;
              break;
          end;
          if (iscell(pspec{3})),
              p2 = p{j};

              for k = 1:length(p2),
                  sz2(1:ndims(p2{k}),k) = size(p2{k})';
              end;
              [match,l,szlist] = matchparamtype(p2,sz2, {pspec{3}}, szlist);
              if (~match),
                  break;
              else
                  done = 1;
                  pspec{end+1} = j;
              end;
          else,
              matchsz = pspec{3};
          end;
         otherwise,
          if (~isa(p{j},pspec{2})),
              match = 0;
              break;
          end;
          matchsz = pspec{3};
        end;

        if (~done),
            matchsz = shiftdim(matchsz);
            if (length(matchsz) == 1),
                p{j} = shiftdim(p{j});
                if (size(p{j},2) ~= 1),
                    match = 0;
                    break
                end;
                sz(:,j) = size(p{j})';
                matchsz(2) = -1;
            end;

            if (ndims(p{j}) ~= length(matchsz)),
                match = 0;
                break;
            else
                if (max(matchsz) > length(szlist)),
                    szlist(max(matchsz),1) = 0;
                end;

                k = find(matchsz > 0);
                checksz(k,1) = szlist(matchsz(k));
                checksz(matchsz < 0) = -matchsz(matchsz < 0);
                checksz(checksz == 0) = sz(checksz == 0,j);

                if (any(checksz ~= sz(:,j))),
                    match = 0;
                    break;
                end;
                szlist(matchsz(k)) = sz(k,j);

                pspec{end+1} = j;
            end;
        end;
        r1{j} = pspec;
    end;

    len = len+1;
end;

