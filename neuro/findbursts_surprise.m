function [burstctr, burstind, burstSurprise, surprise] = findbursts_surprise(t,spike, opt)

spikerate = zeros(size(t));
spikerate(2:end-1) = 1./((t(3:end) - t(1:end-2))/2);

avgspikerate = length(t) / (t(end) - t(1));
%find the shortest length of time in which getting 2 or fewer spikes would
%be unlikely (P < 0.25)
div = 1;
n = Inf;
while (n > 2),
    n = poissinv(0.25,avgspikerate/div);
    div = div+1;
end;
burstEndTime = 1/(div-1);

surpriseThresh = -log(0.05);
maxNotSurprising = 3;

N = length(t);

surprise = zeros(size(t));
b1 = 1;
ind = 1;
while (b1 <= N-3),
    surprise(b1+2) = estSurprise(3,avgspikerate * ...
                                 (t(b1+2)-t(b1)));
    if ((t(b1+2)-t(b1) <= burstEndTime) && ...
        (surprise(b1+2) >= surpriseThresh)),
        nNotSurprising = 0;
        b2 = b1+3;
        while ((b2 < N) && ...
               (nNotSurprising < maxNotSurprising) && ...
               (surprise(b2-1) >= surpriseThresh) && ...
               (t(b2)-t(b2-1) <= burstEndTime)),
            surprise(b2) = estSurprise(b2-b1+1,avgspikerate * ...
                                       (t(b2)-t(b1)));
            if (surprise(b2) < surprise(b2-1)),
                nNotSurprising = nNotSurprising + 1;
            else
                %reset the number of not surprising spikes if the
                %surprise begins increasing again
                nNotSurprising = 0;
            end;
            b2 = b2+1;
        end;

        if ((surprise(b2-1) >= surpriseThresh) && ...
            (b2-b1 >= 5)),
            up(ind) = b1;
            down(ind) = b2-1;
            burstSurprise(ind) = surprise(b2-1);
            ind = ind+1;
        end;

        b1 = b2;
    else
        b1 = b1+1;
    end;
end;

%get the burst center positions
burstctr = zeros(size(up));
for i = 1:length(up),
    burstctr(i) = mean(t(up(i):down(i)));
end;

%now eliminate multiple bursts within the minimum cycle period
if (~isempty(opt.cycledur)),
    good = true(size(up));
    i = 1;
    while (i <= length(up)-1),
        withinper = find(good(i+1:end) & (burstctr(i+1:end) - burstctr(i) < opt.cycledur));
        if (isempty(withinper)),
            i = i+1;
        else
            %add on the first one that's further away than the
            %minimum cycle period
            next = first(good(i+1:end) & (burstctr(i+1:end) - burstctr(i) >= opt.cycledur));
            withinper = [withinper next];
            
            nspike = down(i+withinper) - up(i+withinper);
            [m,ind] = max(nspike);
            
            %eliminate those within a period that have fewer spikes
            %(but hang on to those that have the same number)
            good(i+withinper(nspike < m)) = false;

            i = i+ind;
        end;
    end;
    
    up = up(good);
    down = down(good);
    burstctr = burstctr(good);
    burstSurprise = burstSurprise(good);
end;

burstctr = burstctr';
burstind = [up; down]';
burstSurprise = burstSurprise';

        
