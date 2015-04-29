function [X,coef] = whiten(Y, varargin)
% function [X,coef] = whiten(Y, varargin)

opt.variancefraction = 0.95;

opt = parsevarargin(opt,varargin, 2);

[coef,X,latent] = princomp(Y);

residvar = cumsum(latent) / sum(latent);

ind = last(residvar <= opt.variancefraction);
if (ind < length(residvar))
    ind = ind+1;
end;
good = false(size(latent));
good(1:ind) = true;

X = X(:,good);

Xstd = std(X);

X = bsxfun(@rdivide,X,Xstd);

coef = coef(:,good);
coef = bsxfun(@times,coef,Xstd);

