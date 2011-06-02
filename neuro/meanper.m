function [m,ind] = meanper(per)
% Calculates the mean period, trying to take into account skipped bursts
% See correctBurstFreq.m

if (size(per,1) > size(per,2)),
    per = per';
end;

edges1 = linspace(min(per(:)),max(per(:)),11);
edges1(end) = Inf;
[n,bins] = histc(per,edges1);

[q,ind] = max(n);
m = mean(per(bins == ind));

edges2 = m * [0 0.5 1.5 1.9 2.1 Inf];
[n,bins] = histc(per,edges2);

m = mean([per(bins == 2) 0.5*per(bins == 4)]);
ind = find(bins == 2);


