function C = xcorr1(win, search)
% function C = xcorr1(win, search)
% Partially normalized cross-correlation.

if (length(win) == length(search)),
	w2 = ceil(length(win)/2);
	
	search = shiftdim(search);
	search = [repmat(NaN,[w2 1]); search; repmat(NaN,[w2 1])];
end;

win = shiftdim(win);

winmean = nanmean(win);
win = win - winmean;
winvar = nansum(win.^2);

search = shiftdim(search);
searchmean = nanmean(search);
search = search - searchmean;
searchvar = nansum(search.^2);

S = im2col(search, [length(win) 1], 'sliding');

C = nansum(repmat(win, [1 size(S,2)]) .* S)  ./  sqrt(winvar .* searchvar);

