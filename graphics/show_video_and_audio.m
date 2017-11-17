function [t, img, ts, sound] = show_video_and_audio(filename, varargin)

opt.starttime = 0;
opt = parsevarargin(opt, varargin, 2);

[sound, Fs] = audioread(filename);

ts = (0:length(sound)-1)/Fs;

vid = VideoReader(filename);

clf;
h(1) = axes('Position',[0 0.2 1 0.9]);

h(2) = axes('Position',[0 0.1 1 0.1]);
plot(h(2), ts, sound, 'k-');
axis(h(2), 'tight');

yl = get(h(2), 'YLim');
htime = addplot(h(2), [0; 0], yl, 'r--');

set(h(2), 'Box','on');

if opt.starttime > 0
    vid.CurrentTime = opt.starttime;
end

fr = readFrame(vid);
himg = image(h(1), fr);
axis(h(1), 'tight','equal','ij','off');

fprintf('Click point');
[x,y] = ginputb(1);
x = round(x);
y = round(y);

dur = vid.Duration;
nframes = ceil(dur * vid.FrameRate);

img = zeros(nframes,1);
t = zeros(nframes,1);

img(1) = sum(fr(y,x, :));
t(1) = vid.CurrentTime;

k = 2;
while hasFrame(vid)
    tvid = vid.CurrentTime;
    set(htime, 'XData', [tvid; tvid]);
    fr = readFrame(vid);
    
    set(himg, 'CData',fr);
    
    img(k) = sum(fr(y,x, :));
    t(k) = vid.CurrentTime;
    drawnow;
    
    k = k+1;
end

figure;
plot(t, img, ts,sound);


    
    
    