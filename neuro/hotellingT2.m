function [P,T2,F] = hotellingT2(X,varargin)

opt.alpha = 0.95;
opt.mu = [];

opt = parsevarargin(opt,varargin,2);

[n p] = size(X);

if (isempty(opt.mu))
    opt.mu = zeros(1,p);
end;

m = mean(X);
W = cov(X);

mmu = m - opt.mu;


%%% CONTINUE - check for appropriate size here
T2 = n * mmu * ( W \ mmu' );

F = (n-p)/((n-1)*p) * T2;

P = 1-fcdf(F,p,n-p);

