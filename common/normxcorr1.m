function C = normxcorr1(win, search)

win = shiftdim(win);

winmean = mean(win);
win = win - winmean;
winvar = sum(win.^2);

search = shiftdim(search);
S = im2col(search, [length(win) 1], 'sliding');

Smean = mean(S);
S = S - repmat(Smean, [size(S,1) 1]);
Svar = sum(S.^2);

C = sum(repmat(win, [1 size(S,2)]) .* S)  ./  sqrt(winvar .* Svar);

