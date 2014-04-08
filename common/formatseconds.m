function str = formatseconds(sec, varargin)

assert(sec >= 0);

opt.secformat = '%3$.0fsec';
opt.minformat = '%2$dm:%3$02.0fs';
opt.hourformat = '%1$dh:%2$02dm';
opt = parsevarargin(opt,varargin, 2);

hr = floor(sec/3600);
mn = floor((sec-hr*3600)/60);
rsec = mod(sec, 60);

if (sec < 60)
    str = sprintf(opt.secformat,hr,mn,rsec);
elseif (sec < 3600)
    str = sprintf(opt.minformat,hr,mn,rsec);
else 
    str = sprintf(opt.hourformat,hr,mn,rsec);
end
