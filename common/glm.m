function P = glm(y,factor,model,iscateg,facnames,varnames)
% function P = glm(y,factor,model,[iscateg,facnames,varnames])
% General linear model.  Performs AN(C)OVA/MAN(C)OVA calculations.
% This is a general regression model to determine the effect on
% y of the independent variables in factor.  They may be categorical
% (for ANOVA) or continuous (for regression/ANCOVA) or some
% combination of the two.  Prints an ANOVA table as output.  For
% multivariate situations, it calculates Wilk's lambda for the
% multivariate case, then standard F tests for each univariate case.
% Allows combinations of random and fixed effects models.
%
% Parameters:
% 	y - [N,nvar] matrix of dependent variable(s), where N is the 
%		number of observations and nvar is the number of variables.  
%		Missing	values can be specified by NaN and are removed.
%	factor - [N,nfac] matrix of the factors to test (independent 
%		variable(s)), where nfac is the number of factors.  Each factor
%		can be categorical or continuous.  glm attempts to detect the
%		difference automatically by determining the number of discrete
%		values in each factor.  The parameter iscateg can be used to
%		force a factor to be categorical or continuous.
%	model - String specifying the model to test for significance.  Each
%		successive factor is identified by a letter, starting from 'A'.
%		Interactions are indicated with '*' and a constant can be added.
%		Example: 'constant + A + B + C + A*B + A*C'
%		The error term for testing each effect can be specified after
%		its effect with a '/'.  Note that any specified error term has
%		to be included as an effect, also.
%		Example: 'constant + A + B/A*B + C/A*C + A*B + A*C' would be an
%		appropriate model if A was random and B and C were fixed.
%		Names for the factors can be specified in the parameter facnames
%		and used in the model specification.
%	iscateg - (optional) [nfac,1] vector specifying whether each factor 
%		is categorical or not.  A 1 means the factor is categorical; 0 means
%		continuous; and NaN lets glm figure it out automatically.
%	facnames - (optional) [nfac,1] cell string vector specifying the names
%		of the factors.  These are used in the ANOVA table output and to
%		specify model, above.  If facnames is passed, the names must be used
%		in model.
%	varnames - (optional) [nvar,1] cell string vector specifying the names
%		of the variables.  These are only used in the output table.
%
% Returns P values for the multivariate tests and each univariate one for each
% effect.
%
% To do:
%	- Return a structure appropriate for pairwise testing with multcompare.
%	- Handle nested effects and split plots
%	- Return least square values for each effect
%	- Automatically calculate appropriate error terms for mixed models.

tol = 1e-10;

if (nargin < 6)
	varnames = [];
	if (nargin < 5)
		facnames = [];
		if (nargin < 4),
			iscateg = [];
			if (nargin < 3)
				error('Too few arguments');
			end;
		end;
	end;
end;

nonan = find(all(isfinite(y),2));
if (length(nonan) ~= size(y,1))
	fprintf(1, '%d cases removed due to missing values\n', size(y,1) - length(nonan));
end;

y = y(nonan,:);

% total number of data points
N = size(y,1);
nvar = size(y,2);
if (isempty(varnames)),
	varnames = cellstr(num2str((1:nvar)'));
end;

% number of factors
nfac = size(factor,2);
if (isempty(facnames)),
	facnames = cellstr(char((1:nfac)' + double('A') - 1));
end;

% categorical variables
if (isempty(iscateg))
	iscateg = repmat(NaN,[1 nfac]);
end;

if (~iscell(factor)),
	% determine the number of levels of each factor
	fac = double(factor);
	fac = fac(nonan,:);
	
	sfac = sort(fac);
	dfac = diff(sfac) ~= 0;
	nlev = sum(dfac) + 1;
	
	levelnames = cell(1,nfac);
	
	% assign numbers to each of the factor levels
	q = find(dfac == 1);
	k = 1;
	for i = 1:nfac,	
		% attempt to determine if the factor is categorical
		if (isnan(iscateg(i)) & isnumeric(factor)),
			if (nlev(i) < 0.1*N)
				iscateg = 1;
			else
				iscateg = 0;
			end;
		end;
		
		if (iscateg(i) == 1),
			val = 1;	
			name1 = cell(nlev(i),1);
			for j = 1:nlev(i)-1,
				p = find(fac(:,i) == sfac(q(k)));
				fac(p,i) = val;
				if (isnumeric(factor))
					name1{val} = num2str(factor(p(1),i));
				else
					name1{val} = factor(p(1),i);
				end;
				
				k = k+1;
				val = val+1;
			end;
			p = find(fac(:,i) == sfac(end,i));
			fac(p,i) = val;	
			if (isnumeric(factor))
				name1{val} = num2str(factor(p(1),i));
			else
				name1{val} = factor(p(1),i);
			end;
			
			levelnames{i} = name1;
		else
			nlev(i) = 1;
		end;
	end;
else
	fac = zeros(N,nfac);
	
	for i = 1:nfac,
		if (iscell(factor{i}))
			fac1 = vertcat(factor{i,:});
		else
			fac1 = factor{i};
		end;
		fac1 = fac1(nonan,:);
		f1 = double(fac1);

		% attempt to determine if the factor is categorical
		if (isnan(iscateg(i))),
			if (isnumeric(fac1))
				sfac = sort(fac1);
				dfac = diff(sfac) ~= 0;
				nlev1 = sum(dfac) + 1;
				
				% if the number of levels is less than 10% of the number of 
				% data points, we'll call it categorical
				if (nlev1 < 0.1*N)		
					iscateg(i) = 1;
				else
					iscateg(i) = 0;
				end;
			else
				iscateg(i) = 1;
			end;
		end;
		
		if (iscateg(i) == 1),
			name1 = {};
			val = 1;	
			for j = 1:N,
				if (~any(isnan(f1(j,:)))),
					p = find(all(f1 == repmat(f1(j,:),[N 1]),2));
					
					if (isnumeric(fac1))
						name1{val} = num2str(fac1(p(1),:));
					else
						name1{val} = fac1(p(1),:);
					end;
					fac(p,i) = val;
					f1(p,:) = NaN;
					val = val+1;
				end;
			end;
		
			nlev(i) = val-1;
			levelnames{i} = name1;
		else
			nlev(i) = 1;
			fac(:,i) = fac1;
		end;
	end;
end;

fprintf(1,'Factors\n');
for i = 1:nfac,
	if (iscateg(i))
		levstr = sprintf('%s, ',levelnames{i}{:});
		levstr = levstr(1:end-2);
		
		fprintf(1, '%-10s %d levels (%s)\n', facnames{i}, nlev(i), levstr);
	else
		fprintf(1, '%-10s continuous variable\n', facnames{i});
	end;
end;
fprintf('\n');

% construct the model, if we got text
if (ischar(model)),
	modelstr = lower(model);
	modstr = modelstr;
	
	strrep(modstr,'* ','*');
	strrep(modstr,' *','*');
	
	i = 1;
	model = zeros(1,nfac);
	modlev = zeros(1,nfac+1);
	while (~isempty(modstr)),
		[tok,modstr] = strtok(modstr);
		if (tok == '+'),
			continue;
		end;
		
		[eff,err] = strtok(tok,'/');
		
		model(i,:) = zeros(1,nfac);
		if (strcmp(eff,'constant'))
			modlev(i,end) = 1;
			df(i) = 1;
		else
			k = [];
			for j = 1:nfac,
				if (~isempty(findstr(lower(facnames{j}), eff)))
					model(i,j) = 1;
					k = [k j];
				end;
			end;
			modlev(i,[k end]) = cumprod([1 nlev(k)]);
			
			nlev1 = nlev(k);
			nlev1(nlev1 > 1) = nlev1(nlev1 > 1) - 1;
			df(i) = prod(nlev1);
		end;
		
		i = i+1;
	end;

	nmodel = size(model,1);
	nlevel = sum(modlev(:,end));

	% process error terms
	modstr = modelstr;
	i = 1;
	while (~isempty(modstr))
		[tok,modstr] = strtok(modstr);
		if (tok == '+'),
			continue;
		end;
		
		[eff,err] = strtok(tok,'/');
		if (isempty(err)),
			moderr(i) = 0;
		elseif (strcmp(err,'constant')),
			k = find(all(model == zeros(1,nfac),2));
			
			if (isempty(k))
				error(sprintf('Unknown effect %s for error\n', err));
			end;
			moderr(i) = k;
		else
			bin = zeros(1,nfac);
			for j = 1:nfac,
				if (~isempty(findstr(lower(facnames{j}), err)))
					bin(j) = 1;
				end;
			end;
			k = find(all(model == repmat(bin,[nmodel 1]),2));
			
			if (isempty(k))
				error(sprintf('Unknown effect %s for error\n', err));
			end;
			moderr(i) = k;
		end;
		i = i+1;
	end;
end;

% make dummy variables for each model term
dummy = cell(1,nmodel);
% first do the main effects
k = find(sum(model,2) <= 1);
for i = 1:length(k),
	ifac = find(model(k(i),:) == 1);		% factor number
	
	if (isempty(ifac))
		dummy{k(i)} = ones(N,1);
	else
		dummy{k(i)} = makedummy(fac(:,ifac),nlev(ifac));
	end;
end;

% now do interaction terms
k = find(sum(model,2) > 1);
for i = 1:length(k),
	ifac = find(model(k(i),:) == 1);		% factor numbers
	
	% find each main effect represented in the interaction and
	% cross them
	dum = [];
	for j = 1:length(ifac),
		dum1 = makedummy(fac(:,ifac(j)),nlev(ifac(j)));		
		dum = crossdummy(dum, dum1);
	end;
	
	dummy{k(i)} = dum;
end;

% construct design matrix
% construct model term correspondence matrix modterm
% modterm matrix tells us which column of the design matrix refers to which term
% in the model.  When we do the QR decomp we end up with a reduced design matrix.
% this allows us to selectively remove terms from the reduced matrix
X = zeros(N,sum(df));
modterm = zeros(1,sum(df));
k = 1;
for j = 1:nmodel,
	X(:,k:k+df(j)-1) = dummy{j};
	modterm(k:k+df(j)-1) = j;
	k = k+df(j);
end;

if (rank(X) ~= sum(df))
	error('Oops!');
end;

% solve the system y = X b for the least squares minimized value b

% Tolerance for computing rank from diag(R) after QR decomposition
ncols = size(X,2);
tol = 100 * eps * ncols;

% Do qr decomposition on current design matrix
[q,r,e] = qr(X,0);
if (size(r,2) > 1),
   d = abs(diag(r));
else
   d = abs(r(1,1));
end;
p = sum(d > tol);

% Get information for fit of y on this matrix
yy = (y' * q(:,1:p))';
SSmod = yy'*yy;
SStot = y'*y;

% Return reduced design matrix.
n = size(y,1);
m = min(n, ncols);
X = r(1:p,:);

% Rotate y into the column space of X and keep only that portion.
y = yy(1:min(n,p),:);

% Re-arrange terms to correspond to new design matrix
modterm = modterm(e);

SSerr = SStot - SSmod;

% degrees of freedom
dfErr = N - sum(df);
dfMod = sum(df);

% calculate coeficients for later pairwise comparisons
coeffs = X\y;

% test each effect individually

if (all(moderr == moderr(1)))
	sameerr = 1;
else
	sameerr = 0;
end;

SS = zeros(nmodel,nvar,nvar);
for i = 1:nmodel,
	% construct design matrix without the term we're testing
	k = find(modterm ~= i);
	Xr = X(:,k);
	
	% solve the model
	[q,r,e] = qr(Xr,0);
	if (size(r,2) > 1)
       d = abs(diag(r));
	else
       d = abs(r(1,1));
	end
	p = sum(d > tol);
	
	yy = (y' * q(:,1:p))';
	SS(i,:,:) = SSmod - yy(1:p,:)'*yy(1:p,:);
end;

if (nvar > 1),
	% first print the multivariate statistics
	% Wilks lambda
	fprintf(1,'Wilks lambda\n');
	if (sameerr),
		fprintf(1,'%-20s %10s %10s %10s %10s %10s\n','Effect','V','df1','df2','F','P');
		fprintf(1,'%s\n',repmat('-',[1 75]));
	else
		fprintf(1,'%-20s %-15s %10s %10s %10s %10s %10s\n','Effect','Error','V','df1','df2','F','P');
		fprintf(1,'%s\n',repmat('-',[1 91]));
	end;
	
	lambda = det(SSerr)/det(SSerr + SSmod);
	s = sqrt((nvar^2*dfMod^2 - 4)/(nvar^2 + dfMod^2 - 5));
	y = lambda^(1/s);
	df1 = nvar*dfMod;
	df2 = s*(dfErr - 0.5*(nvar - dfMod + 1)) - 0.5*(nvar*dfMod - 2);
	F = (1-y)/y * df2/df1;
	P(1,1) = 1 - fcdf(F,df1,df2);
	
	stat.lambda(1) = lambda;
	stat.df(1,1) = df1;
	stat.dfErr(1,1) = df2;

	if (sameerr),
		fprintf(1,'%-20s %10.3f %10d %10d %10.3f %10.3f\n','Model',lambda,df1,df2,F,P(1,1));
	else
		fprintf(1,'%-20s %-15s %10.3f %10d %10d %10.3f %10.3f\n','Model','',lambda,df1,df2,F,P(1,1));
	end;
	
	for i = 1:nmodel,
		if (moderr(i) == 0)
			err1 = SSerr;
			dfErr1 = dfErr;
			errName = '';
		else
			err1 = squeeze(SS(moderr(i),:,:));
			dfErr1 = df(moderr(i));
			errName = sprintf('%s', makename(facnames,model(moderr(i),:)));
		end;
	
		lambda = det(err1)/det(err1 + squeeze(SS(i,:,:)));
		if (df(i) > 1),
			s = sqrt((nvar^2*df(i)^2 - 4)/(nvar^2 + df(i)^2 - 5));
		else
			s = 1;
		end;
		y = lambda^(1/s);
		df1 = nvar*df(i);
		df2 = s*(dfErr1 - 0.5*(nvar - df(i) + 1)) - 0.5*(nvar*df(i) - 2);
		F = (1-y)/y * df2/df1;
		P(i+1,1) = 1 - fcdf(F,df1,df2);
	
		stat.lambda(i+1) = lambda;
		stat.F(i+1,1) = F;
		stat.df(i+1,1) = df1;
		stat.dfErr(i+1,1) = df2;
		
		name = makename(facnames,model(i,:));
		if (sameerr)
			fprintf(1,'%-20s %10.3f %10d %10d %10.3f %10.3f\n',name,lambda,df1,df2,F,P(i+1,1));
		else
			fprintf(1,'%-20s %-15s %10.3f %10d %10d %10.3f %10.3f\n',name,errName,lambda,df1,df2,F,P(i+1,1));
		end;
	end;
end;

% univariate tests

for v = 1:nvar,
	fprintf(1,'\n******* Variable %s\n\n',varnames{v});
	
	% print the table
	fprintf(1,'%-20s %10s %10s %10s %10s %10s\n','Effect','SS','df','MS','F','P');
	fprintf(1,'%s\n',repmat('-',[1 75]));
	
	MSmod = SSmod(v,v)/dfMod;
	MSerr = SSerr(v,v)/dfErr;
	F = MSmod/MSerr;
	P(1,v+1) = 1 - fcdf(F,dfMod,dfErr);
	fprintf(1,'%-20s %10.3f %10d %10.3f %10.3f %10.3f\n','Model',SSmod(v,v),dfMod,MSmod,F,P(1,v+1));
	
	for i = 1:nmodel,
		if (moderr(i) == 0)
			err1 = SSerr(v,v);
			dfErr1 = dfErr;
			errName = 'Error (Model)';
		else
			err1 = SS(moderr(i),v,v);
			dfErr1 = df(moderr(i));
			errName = sprintf('Error (%s)', makename(facnames,model(moderr(i),:)));
		end;
		
		name = makename(facnames,model(i,:));
		F = (SS(i,v,v)/df(i))/(err1/dfErr1);
		P(i+1,v+1) = 1 - fcdf(F,df(i),dfErr1);
		
		fprintf(1,'%-20s %10.3f %10d %10.3f %10.3f %10.3f\n',name,SS(i,v,v),df(i),SS(i,v,v)/df(i),F,P(i+1,v+1));
		if (~sameerr),
			fprintf(1,'  %-18s %10.3f %10d %10.3f\n',errName,err1,dfErr1,err1/dfErr1);
		end;

	end;
	
	if (sameerr),
		fprintf(1,'%-20s %10.3f %10d %10.3f\n\n','Error',SSerr(v,v),dfErr,MSerr);
	end;
	
end;

function dum = makedummy(fac, nlev)

if (nlev == 1),		% continuous variable
	dum = fac;
else
	dum = zeros(size(fac,1),nlev-1);
	p = find(fac == 1);
	dum(p,:) = -1;
	
	for j = 2:nlev,
		p = find(fac == j);
		dum(p,j-1) = 1;
	end;
end;

function dum = crossdummy(a,b)

if (isempty(a)), 
	dum = b; 
else
	na = size(a,2);
	nb = size(b,2);
	acols = repmat((1:na), 1, nb);
	bcols = reshape(repmat((1:nb), na, 1), 1, na*nb);
	dum = a(:,acols) .* b(:,bcols);
end;

function name = makename(facnames, model)

name = '';

k = find(model ~= 0);
if (isempty(k))
	name = 'Constant';
else,
	for i = 1:length(k)-1,
		name = [name facnames{k(i)} '*'];
	end;
	name = [name facnames{k(end)}];
end;
