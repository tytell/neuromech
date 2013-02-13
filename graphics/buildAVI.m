function buildAVI(aviname, im, varargin)
% function buildAVI(aviname, im, ...)
% Produces an avi from the 3D image sequence im, where the images are
% sequenced along the third dimension.
% You can provide extra parameters to be passed as options to avifile.
% For example, buildAVI('test.avi',im,'Compression','MSVC','Quality',10).
% Default compression is none.
%
% See avifile for more info.

fig = figure;
set(fig,'DoubleBuffer','on');
set(gca, 'NextPlot', 'replace');

options = varargin;
if (isempty(options)),
	options = {'Compression','none'};
end;

mov = avifile(aviname, options{:});

for i = 1:size(im,3),
	imshow(im(:,:,i));

    F = getframe(gca, [2 2 size(im,2) size(im,1)]);
    mov = addframe(mov,F);
end
mov = close(mov);
delete(fig);
