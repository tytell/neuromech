function setquivercol(h, col)
% function setquivercol(h, col)
%   or     setquivercol(col)
%
% Changes the color of the quiverc plot in h, or finds all plots in the
% current figure and changes their color.  col can be a string or a 3
% element RGB.


if (nargin == 1),
    col = h;
    h = findobj(gcf, 'Tag','quiverc');

    if (isempty(h)),
        error('No quiverc found.');
    end;
end;

if (ischar(col)),
    switch lower(col),
     case {'y','yellow'},
      col = [1 1 0];
     case {'m','magenta'},
      col = [1 0 1];
     case {'c','cyan'},
      col = [0 1 1];
     case {'r','red'},
      col = [1 0 0];
     case {'g','green'},
      col = [0 1 0];
     case {'b','blue'},
      col = [0 0 1];
     case {'w','white'},
      col = [1 1 1];
     case {'k','black'},
      col = [0 0 0];
    end;
end;

for i = 1:length(h),
    switch get(h(i),'Type'),
     case 'line',
      set(h(i),'Color',col);
     case 'patch',
      if (isnumeric(get(h(i),'EdgeColor'))),
          set(h(i),'EdgeColor',col);
      end;
      fcol = get(h(i),'FaceColor');
      if (isnumeric(fcol)),
          set(h(i),'FaceColor',col);
      elseif (~strcmp(fcol,'none')),
          cd = get(h(i),'CData');
          if (prod(size(cd)) == 3),
              set(h(i),'CData',col);
          else
              error(['Quiverc has a specific color per face. ' ...
                     'Cannot override.']);
          end;
      end;
    end;
end;


      
