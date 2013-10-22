function varargout = crosstabG(E,O, varargin)

opt.computemarginal = true;
opt.computeindividual = true;
opt.rownames = {};
opt.colnames = {};
opt.showtable = false;
opt.adjustexpected = false;

sz = size(E);
if ((nargin == 1) || ~isnumeric(O) || any(size(O) ~= sz))
    if (nargin >= 2)
        args = {O varargin{:}};
    else
        args = {};
    end
    O = E;
    E = [];
    p = 2;
else
    p = 3;
    args = varargin;
end

if nargout == 0
    opt.showtable = true;
end

opt = parsevarargin(opt, args, p);

istranspose = false;
nd = ndims(E);
if ((nd == 2) && any(sz == 1))
    if (isempty(E))
        error('1D crosstab does not make sense without an expected frequency');
    end
    nd = 1;
    if (size(O,1) == 1)
        O = O';
        E = E';
        istranspose = true;
    end
    sz = length(O);
end
    
if (nd > 2)
    error('Cannot handle higher than 2D tables yet');
end

S = cell(nd,1);
for i = 1:nd
    S{i} = sum(O,i);
end
N = sum(O(:));
Ne = sum(E(:));

if (nd == 2) && isempty(E)
    %compute expected frequencies
    E = S{2} * S{1} / N;
    
    %doesn't make sense to do the one way if the expected frequencies are
    %not provided
    opt.computemarginal = false;
elseif (abs(Ne - N) > 1)
    if opt.adjustexpected
        E = E/Ne * N;
    else
        warning('crosstabG:badE','Total number expected does not match total observed. G stat may be weird.');
    end
end

if all(sz == 2)
    %Yates correction (Zar, p.505)
    corr = sign(E-O) * 0.5;
    O = O + corr;
end

G = 2 * sum( sum( O .* log(O ./ E)));
df = prod(sz - 1);
P = 1 - chi2cdf(G,df);

Gstat.E = E;
Gstat.O = O;
Gstat.G = G;
Gstat.df = df;
Gstat.P = P;

if nd > 1
    if opt.computemarginal
        [Gstat.Gcol,Gstat.dfcol,Gstat.Pcol] = crosstabG(sum(E,1),sum(O,1), ...
            'computemarginal',false, 'computeindividual',false, ...
            'adjustexpected',opt.adjustexpected);
        [Gstat.Grow,Gstat.dfrow,Gstat.Prow] = crosstabG(sum(E,2),sum(O,2), ...
            'computemarginal',false, 'computeindividual',false, ...
            'adjustexpected',opt.adjustexpected);
    end
    if opt.computeindividual
        Gstat.Gcolindiv = zeros(1,sz(2));
        Gstat.dfcolindiv = zeros(1,sz(2));
        Gstat.Pcolindiv = zeros(1,sz(2));
        for j = 1:sz(2)
            good = ~isnan(E(:,j)) & ~isnan(O(:,j));
            if (sum(good) > 1)
                [Gstat.Gcolindiv(j),Gstat.dfcolindiv(j),Gstat.Pcolindiv(j)] = crosstabG(E(good,j),O(good,j), ...
                    'adjustexpected',opt.adjustexpected);
            else
                Gstat.Gcolindiv(j) = NaN;
                Gstat.dfcolindiv(j) = 0;
                Gstat.Pcolindiv(j) = NaN;
            end
        end
        
        Gstat.Growindiv = zeros(sz(1),1);
        Gstat.dfrowindiv = zeros(sz(1),1);
        Gstat.Prowindiv = zeros(sz(1),1);
        for i = 1:sz(1)
            good = ~isnan(E(i,:)) & ~isnan(O(i,:));
            if (sum(good) > 1)
                [Gstat.Growindiv(i),Gstat.dfrowindiv(i),Gstat.Prowindiv(i)] = crosstabG(E(i,good),O(i,good), ...
                    'adjustexpected',opt.adjustexpected);
            else
                Gstat.Growindiv(i) = NaN;
                Gstat.dfrowindiv(i) = 0;
                Gstat.Prowindiv(i) = NaN;
            end
        end
        
    end
end

if (istranspose)
    Gstat.E = Gstat.E';
    Gstat.O = Gstat.O';
end

if opt.showtable
    fprintf('%20s %6s %6s %6s\n', '', 'G2','df','P');
    fprintf('%20s %6.2f %6d %6.3f\n', 'independence', Gstat.G, Gstat.df, Gstat.P);
    
    if nd > 1
        if opt.computemarginal
            fprintf('%20s %6.2f %6d %6.3f\n', 'rows differ', Gstat.Grow, Gstat.dfrow, Gstat.Prow);
            fprintf('%20s %6.2f %6d %6.3f\n', 'columns differ', Gstat.Gcol, Gstat.dfcol, Gstat.Pcol);
        end
        if opt.computeindividual
            for i = 1:sz(1)
                if isempty(opt.rownames)
                    rowname1 = sprintf('row %d',i);
                else
                    rowname1 = opt.rownames{i};
                end
                fprintf('%20s %6.2f %6d %6.3f\n', rowname1, Gstat.Growindiv(i), Gstat.dfrowindiv(i), Gstat.Prowindiv(i));
            end

            for j = 1:sz(2)
                if isempty(opt.colnames)
                    colname1 = sprintf('col %d',j);
                else
                    colname1 = opt.colnames{j};
                end
                fprintf('%20s %6.2f %6d %6.3f\n', colname1, Gstat.Gcolindiv(j), Gstat.dfcolindiv(j), Gstat.Pcolindiv(j));
            end
        end
    end
end

if nargout == 1
    varargout = {Gstat};
elseif nargout == 3
    varargout = {G, df, P};
end

    

