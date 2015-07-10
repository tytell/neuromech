function [step,direc,enab] = stepper(t,pos, varargin)

opt.outfreq = 100000;   % Hz
opt.stepsperrev = 4000;
opt.enablehi = false;

opt = parsevarargin(opt,varargin, 3);

stepsize = 360/opt.stepsperrev;

t2 = (t(1):1/opt.outfreq:t(end))';
pos2 = interp1(t,pos, t2);

step = false(size(t2));
direc = false(size(t2));
if opt.enablehi
    enab = ones(size(t2));
else
    enab = zeros(size(t2));
end

enab([1 end]) = 1-enab([1 end]);

d = [0; diff(pos2)];
direc = d >= 0;

prevpos = pos(1);
for i = 2:length(pos2)
    if step(i-1)
        continue;
    elseif pos2(i) - prevpos >= stepsize
        step(i) = 1;
        prevpos = prevpos + stepsize;
    elseif pos2(i) - prevpos <= -stepsize
        step(i) = 1;
        prevpos = prevpos - stepsize;
    end
end

        