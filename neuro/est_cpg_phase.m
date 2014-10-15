function [perstart,burstperstart,burstcyclenum] = est_cpg_phase(t,burstt,chanlocation, varargin)

opt.tol = 0.1;
opt.smoothwindow = 9;
opt.outlier = 1.5;
opt.relativeto = 'L1';
opt.basephaselag = 0.01;
opt.excludenear = [];
opt.excludenum = 0;

opt = parsevarargin(opt,varargin,4);

nchan = size(burstt,2);

side = regexp(chanlocation,'[LR]','once','match');
isleft = cellfun(@(x) x == 'L', side);

seg = regexp(chanlocation,'\d+','once','match');
seg = cellfun(@str2double,seg);

relside = regexp(opt.relativeto,'[LR]','once','match');
relseg = regexp(opt.relativeto,'\d+','once','match');
relseg = str2double(relseg);

burstrelphase = (seg - relseg)*opt.basephaselag;
if (relside == 'L')
    burstrelphase(~isleft) = burstrelphase(~isleft) + 0.5;
else
    burstrelphase(isleft) = burstrelphase(isleft) + 0.5;
end

burstper = zeros(size(burstt));
for i = 1:nchan
    good = isfinite(burstt(:,i));
    burstt1 = burstt(good,i);
    per1 = NaN(size(burstt1));
    per1(2:end-1) = (burstt1(3:end) - burstt1(1:end-2))/2;
    burstper(good,i) = per1;
end

good = isfinite(burstt) & isfinite(burstper);
t2 = burstt(good);
freq2 = 1./burstper(good);
if (size(freq2,2) == 1)
    freq2 = freq2';
    t2 = t2';
end

[t2,ord] = sort(t2);
freq2 = freq2(ord);
revord(ord) = 1:length(freq2);

halfwind = floor(opt.smoothwindow/2);

fwind = vec2col(freq2,2*halfwind+1,1);
dropmiddle = [1:halfwind halfwind+2:size(fwind,1)];

fiqr = iqr(fwind(dropmiddle,:));
fmed = median(fwind(dropmiddle,:));

isoutlier = abs(freq2 - fmed) > opt.outlier * fiqr;

freq2(isoutlier) = NaN;

fwind = vec2col(freq2,2*halfwind+1,1);

sysfreq2 = nanmedian(fwind);
sysphase2 = cumtrapz(t2,sysfreq2);

perstart0 = interp1(sysphase2,t2, ceil(min(sysphase2)):floor(max(sysphase2)));

burstcyclenum = NaN(size(burstt));
burstperstart0 = NaN(size(burstt));
for c = 1:nchan
    j = 1;
    for i = 1:size(burstt,1)
        if isfinite(burstt(i,c))
            while (j < length(perstart0)) && (perstart0(j+1) <= burstt(i,c))
                j = j+1;
            end
            burstcyclenum(i,c) = j;
            burstperstart0(i,c) = perstart0(j);
        end
    end
end

burstphase = (burstt - burstperstart0) ./ burstper;
cpgphase = burstphase - repmat(burstrelphase(:)',[size(burstphase,1) 1]);

good = isfinite(burstcyclenum);
cpgphasemn = accumarray(burstcyclenum(good),cpgphase(good),[],@phasemean1);
cpgtmn = accumarray(burstcyclenum(good),burstt(good),[],@nanmean);

perstart = interp1(cpgphasemn + (1:length(cpgphasemn))', cpgtmn, (1:length(cpgphasemn))');

[burstperstart,burstcyclenum] = get_cycle_time(perstart, burstt);

function pmn = phasemean1(x)

C = cos(2*pi*x);
S = sin(2*pi*x);
n = sum(isfinite(x));

pmn = atan2(nanmedian(S),nanmedian(C)) ./ (2*pi);
pmn = mod(pmn,1);





