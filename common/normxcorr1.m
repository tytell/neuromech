function [C,xs] = normxcorr1(win, search)
% function [C,xs] = normxcorr1(win, search)

if (length(win) == length(search)),
	w2 = ceil(length(win)/2);
	
	search = shiftdim(search);
	search = [repmat(NaN,[w2 1]); search; repmat(NaN,[w2 1])];
end;

w2 = (length(search)-length(win))/2;
xs = -ceil(w2):floor(w2);

win = shiftdim(win);

winmean = nanmean(win);
win = win - winmean;
winvar = nansum(win.^2);

search = shiftdim(search);
S = im2col(search, [length(win) 1], 'sliding');

Smean = nanmean(S);
S = S - repmat(Smean, [size(S,1) 1]);
Svar = nansum(S.^2);

C = nansum(repmat(win, [1 size(S,2)]) .* S)  ./  sqrt(winvar .* Svar);

