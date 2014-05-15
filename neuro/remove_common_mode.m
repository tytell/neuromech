function [rem,common,commonmag] = remove_common_mode(sig, varargin)

opt.clip = 10;
opt.average = 'median';

sig(abs(sig) >= opt.clip) = NaN;

switch opt.average
    case 'mean'
        common = nanmean(sig,2);
    case 'median'
        common = nanmedian(sig,2);
    otherwise
        error('Unrecognized average method');
end

common = common ./ sqrt(nanmean(common.^2));

rem = sig;
commonmag = zeros(1,size(sig,2));
for i = 1:size(sig,2)
    commonmag(i) = sqrt(nanmean((sig(:,i).*common).^2));
    
    rem(:,i) = rem(:,i) - common.*commonmag(i);
end


