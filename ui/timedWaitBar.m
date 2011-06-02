function isGood = timedWaitBar(pct, varargin)
% function isGood = timedWaitBar(pct)
%		or isGood = timedWaitBar(pct, title)
%
% You must call timedWaitBar(0,...) first to initialize the timer.  After
% that, it assumes that the percentage passed is proportional to the amount
% of time, and uses that to calculate the time remaining.  To close the wait
% bar, call it with pct = 1.
%
% If the user presses the cancel button, isGood is false.
%
% Example:
%
% timedWaitBar(0,'Calculating...');
% for i = 1:20,
%	% do calculation
%	if (~timedWaitBar(i/20))
%		break;
%	end;
% end;
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell

global timedWaitBarHandle;
global timedWaitBarStart;

global timedWaitBarDebug;
global timedWaitBarElapsed;
global timedWaitBarStep;
global timedWaitBarError;
global timedWaitBarPrev10;

if (ischar(pct)),
    switch lower(pct),
     case 'cancel',
      if (ishandle(timedWaitBarHandle)),
          delete(timedWaitBarHandle);
      end;
      timedWaitBarHandle = -1;
      timedWaitBarStart = -1;
      timedWaitBarStep = 1;
    end;
    
    return;
end;

k = find(cellfun('isclass',varargin,'char'));
for i = k,
    s = varargin{k};
    q = find(s == '\');
    s(q) = '/';
    varargin{k} = s;
end;

isGood = 1;
if (pct == 0),				% initialize
    if (strcmp(varargin{1},'debug')),
        timedWaitBarDebug = 1;
        varargin = varargin(2:end);
    else
        timedWaitBarDebug = 0;
        timedWaitBarError = 0;
    end;

    if (isempty(timedWaitBarHandle) | ~ishandle(timedWaitBarHandle)),
        timedWaitBarHandle = waitbar(0,varargin{:},'CreateCancelBtn', ...
                                     'timedWaitBar(''cancel'')');
    else
        waitbar(0,timedWaitBarHandle,varargin{:});
    end;
    
    set(timedWaitBarHandle, 'Name', '');
    
    timedWaitBarStart = cputime;
    timedWaitBarElapsed = [];
    timedWaitBarStep = 1;
    timedWaitBarError = 0;
    timedWaitBarPrev10 = 0;
else
    if ((pct < 1) & ~ishandle(timedWaitBarHandle)),
        if (nargout == 0),
            error('Timed wait bar uninitialized.');
        else
            isGood = 0;
            return;
        end;
    end;
    
    if (pct >= 1),
        if (ishandle(timedWaitBarHandle)),
            delete(timedWaitBarHandle);
        end;
        timedWaitBarHandle = -1;
        timedWaitBarStart = -1;
    else,
        elapsed = cputime - timedWaitBarStart;
        remain = elapsed/pct * (1-pct);
        
        if (timedWaitBarDebug),
            timedWaitBarElapsed(timedWaitBarStep) = elapsed;
            timedWaitBarStep = timedWaitBarStep + 1;
            
% test every 10 seconds
            if (elapsed-timedWaitBarPrev10 > 10),
                dt = diff(timedWaitBarElapsed)';
                step = (1:timedWaitBarStep-2)';
                
                [b,bint] = regress(dt, [step ones(size(step))]);
                if (sign(bint(1,1)) == sign(bint(1,2))),
                    timedWaitBarError = 1;
                else
                    timedWaitBarError = 0;
                end;
                
                timedWaitBarPrev10 = elapsed;
            end;
        end;
        
        elapStr = fmtTime(elapsed);
        remStr = fmtTime(remain);
        
        if (~timedWaitBarError),
          nameStr = sprintf('%s elapsed, %s remaining.',elapStr,remStr);
        else
          nameStr = sprintf('%s elapsed, %s? remaining.',elapStr,remStr);
        end;
        waitbar(pct,timedWaitBarHandle,varargin{:});
        if (~ishandle(timedWaitBarHandle)),
            isGood = 0;
            return;
        else
            set(timedWaitBarHandle,'Name',nameStr);
        end;
    end;
end;

function str = fmtTime(sec)

min = floor(sec/60);
hour = floor(min/60);
sec = mod(sec,60);

if (hour > 0),
    min = mod(min,60);
    str = sprintf('%dh%02dm',hour,min);
elseif (min > 0),
    str = sprintf('%dm%02ds',min,round(sec));
else
    str = sprintf('%dsec',round(sec));
end;

	
