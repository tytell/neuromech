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
opt.path = 'head';  % or 'COM'
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

if isfield(K,'mxmm')
    npt = size(K.mxmm,1);
elseif isfield(K,'npts')
    npt = K.npts;
end
    
isPeak = zeros([nfr npt]);
for i = 1:npt,
    k = find(isfinite(K.indpeak(i,:)));
    isPeak(K.indpeak(i,k),i) = 1;
end;
tpeak = NaN([size(K.indpeak,2) 1]);
tpeak(k) = t(K.indpeak(end,k));


dt = t(2)-t(1);
assert(isfinite(dt));

if ~isfield(K,'pathx')
    switch opt.path
        case 'head'
            fishpath = get_fish_path(K.t,K.hxmm,K.hymm,K.humms,K.hvmms,K.haxmmss,K.haymmss, ...
                'path','head');
            
            pathcurve = fishpath.pathcurve;
            pathang = fishpath.pathang;
            speed = fishpath.speed;
            accel = fishpath.accel;
                        
        case {'COM','com'}
            fishpath = get_fish_path(K.t,K.mxmm,K.mymm, 'path','COM');
            
            pathcurve = fishpath.pathcurve;
            pathang = fishpath.pathang;
            speed = fishpath.speed;
            accel = fishpath.accel;
    end
else
    pathcurve = K.pathcurve;
    pathang = K.pathang;
    speed = K.speed;
    accel = K.accel;
end
        
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

if isfield(K,'mxmm')
    s = [zeros(1,size(K.mxmm,2)); cumsum(sqrt(diff(K.mxmm).^2 + diff(K.mymm).^2))];
    s = nanmedian(s,2);
else
    s = linspace(0,K.fishlenmm,npt)';
end

tailamp = K.amp(end,:)';
headamp = K.amp(1,:)';
[pkamp,pkamppt] = max(K.amp);
pkamp = pkamp';
pkamppos = s(pkamppt) / s(end);

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
fprintf(fid, '%%%% Length %f\n', K.fishlenmm);
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

if isfield(K,'mxmm')
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
end