function saveKinematicsCSV(outfile,kinfile, varargin)
% function saveKinematicsCSV(outfile,kinfile)
% Converts kinematic data in kinfile (a .mat file) and saves in CSV form to
% outfile
%
% Copyright (c) 2013, Eric Tytell

opt.smoothper = 1;      % in sec
opt.smoothmethod = 'loess';
opt.ampthresh = 0.01;   % fraction of body length
opt.npt = 20;
opt = parsevarargin(opt,varargin, 3);

K = load(kinfile);

if (isfield(K,'frkin')),
    fr = shiftdim(K.frkin);
else
    fr = shiftdim(K.fr);
end;

nfr = length(fr);
nbt = size(K.indpeak,2);

t = shiftdim(K.t);

npt = size(K.mxs,1);
MXY(:,1:2:2*npt) = K.mxmm';
MXY(:,2:2:2*npt) = K.mymm';

H = [K.humms' K.hvmms' K.haxmmss' K.haymmss'];

isPeak = zeros([nfr npt]);
for i = 1:npt,
    k = find(isfinite(K.indpeak(i,:)));
    isPeak(K.indpeak(i,k),i) = 1;
end;
tpeak = NaN([size(K.indpeak,2) 1]);
tpeak(k) = t(K.indpeak(end,k));


dt = t(2)-t(1);
assert(isfinite(dt));

good = isfinite(K.hxmm) & isfinite(K.hymm);
smoothfrac = opt.smoothper / (sum(good)*dt);
pathx = NaN(size(K.hxmm));
pathy = NaN(size(K.hymm));

if (smoothfrac < 1)
    pathx(good) = smooth(K.hxmm(good),smoothfrac,opt.smoothmethod);
    pathy(good) = smooth(K.hymm(good),smoothfrac,opt.smoothmethod);
else
    p = polyfit(t,K.hxmm(good),2);
    pathx(good) = polyval(t(good),p);
    p = polyfit(t,K.hymm(good),2);
    pathy(good) = polyval(t(good),p);
end

pathcurve = curvature(pathx,pathy, 'smooth',1, 'splineindiv');

swimvecx = [diff(pathx) NaN];
swimvecy = [diff(pathy) NaN];
mag = sqrt(swimvecx.^2 + swimvecy.^2);
swimvecx = swimvecx ./ mag;
swimvecy = swimvecy ./ mag;

pathang = unwrap(atan2(swimvecy,swimvecx)+pi) - pi;

speed = K.humms.*swimvecx + K.hvmms.*swimvecy;
accel = K.haxmmss.*swimvecx + K.haymmss.*swimvecy;

speedmn = NaN(nbt,1);
accelmn = NaN(nbt,1);
pathangmn = NaN(nbt,1);
pathcurvemn = NaN(nbt,1);
for j = 1:nbt-1
    if (isfinite(K.indpeak(end,j)) && isfinite(K.indpeak(end,j+1)))
        k = K.indpeak(end,j):K.indpeak(end,j+1);
        
        speedmn(j) = nanmean(speed(k));
        accelmn(j) = nanmean(accel(k));
        pathangmn(j) = angmean(pathang(k));
        pathcurvemn(j) = nanmean(pathcurve(k));
    end
end
        
freq = (1./nanmedian(K.per(end-4:end,:)))';

s = [zeros(1,size(K.mxmm,2)); cumsum(sqrt(diff(K.mxmm).^2 + diff(K.mymm).^2))];
s = nanmedian(s,2);

tailamp = K.amp(end,:)';
headamp = K.amp(1,:)';
[pkamp,pkamppt] = max(K.amp);
pkamp = pkamp';
pkamppos = s(pkamppt) / s(end);

exc = nanmean(K.exc,3)';
if (~isempty(K.wavevel))
    wavevel = K.wavevel';
else
    wavevel = NaN(size(tpeak));
end
if (~isempty(K.wavelen))
    wavelen = nanmedian(K.wavelen)';
else
    wavelen = NaN(size(tpeak));
end

fid = fopen(outfile,'w');
fprintf(fid, '%%%% From data file %s: %s\n', kinfile, date);
fprintf(fid, '%%%% Length %f\n', nanmedian(K.smm(end,:)));
fprintf(fid,'\n');

fprintf(fid,'%%%% General kinematics:\n');
ttl = {'PeakTime(s)','Speed(mm/s)','Accel(mm/s^2)','PathAng(deg)','PathCurve(1/mm)',...
    'BeatFreq(Hz)','HeadAmp(mm)','TailAmp(mm)','PeakAmp(mm)', ...
              'PeakAmpPos(%)','WaveLen(mm)','WaveVel(mm/s)'};
fprintf(fid, '%s,', ttl{:});
fprintf(fid, '\n');

Kin1 = [tpeak speedmn accelmn pathangmn*180/pi pathcurvemn ...
    freq headamp tailamp pkamp pkamppos*100 wavelen wavevel];
tplt = repmat('%10.4f,',[1 size(Kin1,2)]);
tplt = [tplt(1:end-1) '\n'];
fprintf(fid, tplt, Kin1');

fclose(fid);

[pn,fn,ext] = fileparts(outfile);
outfile2 = fullfile(pn,[fn '-points' ext]);

fid = fopen(outfile2,'w');
fprintf(fid, '%%%% From data file %s: %s\n', kinfile, date);
fprintf(fid, '%%%% Length %f\n', nanmedian(K.smm(end,:)));
fprintf(fid,'\n');

fprintf(fid, '%%%% Midline coordinates (mm)\n');
ttl = {'FrameNum','Time(s)','PointNum','X(mm)','Y(mm)'};
fprintf(fid, '%s,', ttl{:});
fprintf(fid, '\n');

[ptnum,frnum] = ndgrid(1:K.npts, K.fr);
t2 = repmat(K.t,[K.npts 1]);
Pts1 = [frnum(:) t2(:) ptnum(:) K.mxmm(:) K.mymm(:)];
tplt = repmat('%10.4f,',[1 size(Pts1,2)]);
tplt = [tplt(1:end-1) '\n'];
fprintf(fid, tplt, Pts1');

fclose(fid);
