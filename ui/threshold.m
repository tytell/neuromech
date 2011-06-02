function t = threshold(I)
% function t = threshold(I)

done = 0;

clf;
imshow(I,'n');

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
      for i = i:size(tmap,1),
         tmap(i,:) = [1,0,0];
      end;
      colormap(tmap);
   end;
end;

colormap(old);

         
  