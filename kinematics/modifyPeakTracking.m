function indpeak = modifyPeakTracking(indpeak)

s = 1:size(indpeak,1);

hand = plot(s,indpeak);
xlabel('Body point');
ylabel('Frame in which peak was found');

set(hand,'ButtonDownFcn',{@peakClick_callback,hand});

fig = gcf;
set(fig,'KeyPressFcn',@keyPress_callback, 'UserData','');

done = 0;
while ~done,
    waitfor(fig,'UserData');

    pk = get(fig,'UserData');

    if (pk == 0),
        done = 1;
    else
        act = input('Action? (d - delete, l - remove left, q - quit) ','s');
        if (~isempty(act)),
            switch lower(act),
             case 'd',
              indpeak(:,pk) = NaN;
             case 'l',
              disp('Click first good point.');
              [pt,q] = ginput(1);
              indpeak(1:floor(pt),pk) = NaN;
             case 'q',
              done = 1;
            end;
        end;

        set(hand(pk),'YData',indpeak(:,pk));
    end;
    set(fig,'UserData','');
end;

k = find(any(isfinite(indpeak)));
indpeak = indpeak(:,k);


% ------------------------------------------------------
function peakClick_callback(obj, eventdata, handles)

set(handles,'LineWidth',0.5,'Marker','none');
set(obj,'LineWidth',2,'Marker','.');

fig = get(get(obj,'Parent'),'Parent');
k = find(handles == obj);
set(fig,'UserData',k);



% ------------------------------------------------------
function keyPress_callback(obj, eventdata)

c = get(obj,'CurrentCharacter');

switch lower(c),
 case {'q',char(13)},                   % q or return quits
  set(obj,'UserData',0);
end;

