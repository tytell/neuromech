function saveKinematicsCSV(outfile,kinfile)
% function saveKinematicsCSV(outfile,kinfile)
% Converts kinematic data in kinfile (a .mat file) and saves in CSV form to
% outfile
%
% Mercurial revision hash: $Revision: dd099ded8eec $ $Date: 2010/08/11 13:58:29 $
% Copyright (c) 2010, Eric Tytell

K = load(kinfile);

if (isfield(K,'frkin')),
    fr = shiftdim(K.frkin);
else
    fr = shiftdim(K.fr);
end;

nfr = length(fr);

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

freq = (1./nanmedian(K.per(end-4:end,:)))';

s = [zeros(1,size(K.mxmm,2)); cumsum(sqrt(diff(K.mxmm).^2 + diff(K.mymm).^2))];
s = nanmedian(s,2);

tailamp = K.amp(end,:)';
headamp = K.amp(1,:)';
[pkamp,pkamppt] = max(K.amp);
pkamp = pkamp';
pkamppos = s(pkamppt) / s(end);

exc = nanmean(K.exc,3)';
wavevel = K.wavevel';
wavelen = nanmedian(K.wavelen)';

fid = fopen(outfile,'w');
fprintf(fid, '%%%% From data file %s: %s\n', kinfile, date);
fprintf(fid,'\n');

fprintf(fid,'General kinematics:\n');
ttl = {'PeakTime(s)','BeatFreq(Hz)','HeadAmp(mm)','TailAmp(mm)','PeakAmp(mm)', ...
              'PeakAmpPos(%)','WaveLen(mm)','WaveVel(mm/s)'};
fprintf(fid, '%s,', ttl{:});
fprintf(fid, '\n');

Kin1 = [tpeak freq headamp tailamp pkamp pkamppos*100 wavelen wavevel];
tplt = repmat('%10.4f,',[1 size(Kin1,2)]);
tplt = [tplt(1:end-1) '\n'];
fprintf(fid, tplt, Kin1');

fprintf(fid,'\n\nAmplitude Envelope (mm):\n');
ttl = cell(1,npt+1);
ttl{1} = 'PeakTime(s)';
for i = 1:npt,
    ttl{i+1} = num2str(i);
end;
ttl{2} = 'Head';
ttl{end} = 'Tail';

fprintf(fid, '%s,', ttl{:});
fprintf(fid, '\n');

Kin1 = [tpeak exc];
tplt = repmat('%10.4f,',[1 size(Kin1,2)]);
tplt = [tplt(1:end-1) '\n'];
fprintf(fid, tplt, Kin1');

fprintf(fid,'\n\nVelocity and Acceleration:\n');

ttl = {'Frame','Time(s)','HeadVelX(mm/s)','HeadVelY(mm/s)','HeadAccX(mm/s)','HeadAccY(mm/s)'};
fprintf(fid, '%s,', ttl{:});
fprintf(fid, '\n');

Kin1 = [fr t H];
tplt = repmat('%10.4f,',[1 size(Kin1,2)]);
tplt = [tplt(1:end-1) '\n'];
fprintf(fid, tplt, Kin1');

fclose(fid);
