function [cyclet,cycleind] = get_cycle_time(cyclestart, t)

[~,ord] = sort([cyclestart(:); t(:)]);
inc = [ones(numel(cyclestart),1); zeros(numel(t),1)];
inc = inc(ord);

revord(ord) = 1:length(inc);

inc(1) = 1;
cycleind = cumsum(inc);
cycleind = cycleind(revord);
cycleind = cycleind(numel(cyclestart)+1:end);
cycleind(cycleind == 0) = 1;
cycleind(cycleind > numel(cyclestart)) = numel(cyclestart);

cyclet = cyclestart(cycleind);

cyclet = reshape(cyclet, size(t));
cycleind = reshape(cycleind, size(t));

