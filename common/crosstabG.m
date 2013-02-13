function varargout = crosstabG(E,O, varargin)

opt.computemarginal = true;
opt.computeindividual = true;

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

if (nd == 2) && isempty(E)
    %compute expected frequencies
    E = S{2} * S{1} / N;
    
    %doesn't make sense to do the one way if the expected frequencies are
    %not provided
    opt.computemarginal = false;
end

if all(sz == 2)
    %Yates correction
    corr = sign(E-O) * 0.5;
    E = E + corr;
end

G = 2 * sum( sum( O .* log(O ./ E)));
df = prod(sz - 1);
P = 1 - chi2cdf(G,df);

Gstat.E = E;
Gstat.O = O;
Gstat.df = df;
Gstat.P = P;

if nd > 1
    if opt.computemarginal
        [Gstat.Gcol,Gstat.dfcol,Gstat.Pcol] = crosstabG(sum(E,1),sum(O,1), ...
            'computemarginal',false, 'computeindividual',false);
        [Gstat.Grow,Gstat.dfrow,Gstat.Prow] = crosstabG(sum(E,2),sum(O,2), ...
            'computemarginal',false, 'computeindividual',false);
    end
    if opt.computeindividual
        Gstat.Gcolindiv = zeros(1,sz(2));
        Gstat.dfcolindiv = zeros(1,sz(2));
        Gstat.Pcolindiv = zeros(1,sz(2));
        for j = 1:sz(2)
            [Gstat.Gcolindiv(j),Gstat.dfcolindiv(j),Gstat.Pcolindiv(j)] = crosstabG(E(:,j),O(:,j));
        end
        
        Gstat.Growindiv = zeros(sz(1),1);
        Gstat.dfrowindiv = zeros(sz(1),1);
        Gstat.Prowindiv = zeros(sz(1),1);
        for i = 1:sz(1)
            [Gstat.Growindiv(i),Gstat.dfrowindiv(i),Gstat.Prowindiv(i)] = crosstabG(E(i,:),O(i,:));
        end
        
    end
end

if (istranspose)
    Gstat.E = Gstat.E';
    Gstat.O = Gstat.O';
end

if nargout == 1
    varargout = {Gstat};
elseif nargout == 3
    varargout = {G, df, P};
end

    

