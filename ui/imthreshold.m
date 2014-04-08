function t = imthreshold(I,varargin)
% function t = imthreshold(I)

opt.direction = 'up';
opt = parsevarargin(opt,varargin,2);

done = 0;

clf;
imshow(I,'InitialMagnification','fit');

old = colormap;
colormap(gray);
map = colormap;

bar = colorbar;
lim = get(bar,'YLim');

while ~done,
    if (waitforbuttonpress == 1)
        done = 1;
    elseif (gca == bar)
        pt = get(gca,'CurrentPoint');
        t = pt(1,2);
        
        tmap = map;
        i = (t-lim(1))/(lim(2)-lim(1));
        i = round(i*size(tmap,1));
        
        switch opt.direction
            case 'up'
                k = i:size(tmap,1);
            case 'down'
                k = 1:i;
        end
        tmap(k,:) = repmat([1,0,0], [length(k) 1]);
        colormap(tmap);
    end;
end;

colormap(old);

         
  