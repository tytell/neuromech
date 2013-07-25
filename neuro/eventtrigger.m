function varargout = eventtrigger(pulseon, tevent, rng, varargin)
% function [A,B,C...] = eventtrigger(pulseon, tevent, rng, A,B,C,...,
% options)
%   Reshapes tevent around pulseon for a certain number of events before
%   and after the pulse.  Optionally also reshapes other variables
%   (A,B,C,...) that have the same size as tevent.  
%
%   Options: 'timerange' - rng specifies a range of absolute time, relative
%                to the pulse
%            'eventrange' - (default) rng species the number of events
%                before and after the pulse
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell <tytell at jhu dot edu>

npulse = length(pulseon);
nchan = size(tevent,1);

i = 1;
while ((i <= length(varargin)) && ~ischar(varargin{i}) && ...
        (ndims(varargin{i}) == ndims(tevent)) && ...
        all(size(varargin{i}) == size(tevent))),
    i = i+1;
end;
extra = varargin(1:i-1);

opt.rangetype = 'eventrange';
opt.method = 'sort';
opt.out = 'matrixout';

opt = parsevarargin(opt,varargin(i:end), ...
    'multival',{'rangetype',{'timerange','eventrange'}, ...
                'method',{'sort','loop'}, ...
                'out',{'cellout','matrixout'}});

switch opt.rangetype,
    case 'eventrange'
        iseventrange = true;
    case 'timerange'
        iseventrange = false;
end;

if (isempty(rng)),
    rng = [-Inf Inf];
    opt.method = 'sort';
end;

if ((ndims(pulseon) > 2) || all(size(pulseon) > 1))
    error('Pulse times must be a vector');
end;
if (size(pulseon,2) == 1),
    %column vector
    pulseon = pulseon';
end;

switch opt.method,
  case 'loop',
    %build up the matrix that defines the reshaped data
    if (iseventrange),
        nevent = rng(2) - rng(1) + 1;
        
        %for an event range, we know how many events will be represented, so
        %we can preallocate the arrays
        event = zeros(nevent,npulse,nchan);
    else
        %for a time window, we don't know, so we have to set the size to zero
        event = zeros(0,npulse,nchan);
    end;

    if (any(flatten(diff(tevent,[],2)) < 0) || any(flatten(diff(pulseon)) < 0)),
        %warning('Works best with monotonically increasing event times.');
        isinc = false;
    else
        isinc = true;
    end;

    for ch = 1:nchan,
        eventind = 1;
        for p = 1:npulse,
            if (iseventrange),
                %find the first event that starts after the pulse came on
                firstevent = first(tevent(ch,:) >= pulseon(p), 2, eventind);
                
                if (~isempty(firstevent)),
                    %and take the range of events
                    rng1 = firstevent + (rng(1):rng(2));

                    event(:,p,ch) = rng1;
                    
                    if (isinc),
                        eventind = firstevent+1;
                    end;
                end;
            else
                %find the range of events that occurred within the time
                %specifed
                a = first(tevent(ch,:) >= pulseon(p)+rng(1), 2, eventind);
                if (isempty(a)),
                    continue;
                end;
                eventind = first(tevent(ch,:) >= pulseon(p), 2, a);
                if (isempty(eventind)),
                    continue;
                end;
                b = first(tevent(ch,:) >= pulseon(p)+rng(2)) - 1;
                if (isempty(b)),
                    continue;
                end;
                
                event(1:(b-a+1),p,ch) = (a:b)';
                
                if (isinc),
                    eventind = eventind + 1;
                else
                    eventind = 1;
                end;
            end;
        end;
    end;

    sz = size(event);
    channel = repmat(shiftdim(1:nchan,-1), [sz(1) npulse 1]);

    %N is the data in tevent, plus any extra
    N = cell(1,length(extra)+1);
    [N{:}] = deal(nans(sz));

    %build indices into the tevent matrix
    good = event ~= 0;
    ind = sub2ind(size(tevent), channel(good), event(good));

    %and do the data rearranging
    N{1}(good) = tevent(ind);
    for i = 2:length(N),
        N{i}(good) = extra{i-1}(ind);
    end;

    varargout = N;

  case 'sort',
    nrng = diff(rng)+1;

    N = cell(1,length(extra)+1);
    if (iseventrange),
        [N{:}] = deal(NaN(nrng,npulse,nchan));
    else
        [N{:}] = deal(cell(npulse,nchan));
    end;
    
    goodpulse = isfinite(pulseon);
    ngoodpulse = sum(goodpulse);
    for ch = 1:nchan,
        goodevent = isfinite(tevent(ch,:));
        nevent = sum(goodevent);

        %sort the times of the pulses in with the times of the events
        [events,ord] = sort(cat(2,pulseon(goodpulse),tevent(ch,goodevent)));
        
        %use the ord matrix to identify which events are pulses and which are not
        ispulse = cat(2,true(1,ngoodpulse),false(1,nevent));
        ispulse = ispulse(ord);
        
        %also make sure we know the order the pulses were in - they might not have been
        %sorted originally, and we want to return them in the same order
        pn1 = 1:length(pulseon);
        pulsenum = cat(2,pn1(goodpulse),zeros(1,nevent));
        pulsenum = pulsenum(ord);
        
        %this gives us the sorted order of the pulses
        pulseord = pulsenum(ispulse);
        
        %and this is the index of the pulse in the interleaved data set
        pulseind = find(ispulse);
        %subtract off the number of interleaved pulses, so that we can reference just
        %the event time, rather than the interleaved data set
        pulseind = pulseind - (0:ngoodpulse-1);
        
        %number of events between each pulse
        eventperpulse = diff([1 pulseind sum(goodevent)+1]);

        if (iseventrange),
            %construct the index for events preceeding each pulse
            indpre = (rng(1):-1)';
            indpre = indpre(:,ones(1,ngoodpulse));
            %don't back up over the previous pulse
            indpre(indpre < -eventperpulse(ones(-rng(1),1),1:end-1)) = NaN;
            indpre = indpre + pulseind(ones(-rng(1),1),:);
            
            %index for events after each pulse
            indpost = (0:rng(2))';
            indpost = indpost(:,ones(1,ngoodpulse));
            indpost(indpost >= eventperpulse(ones(rng(2)+1,1),2:end)) = NaN;
            indpost = indpost + pulseind(ones(rng(2)+1,1),:);
            
            %put the indices together
            ind = NaN(nrng,npulse);
            ind(1:(-rng(1)),pulseord) = indpre;
            ind((-rng(1)+1):end,pulseord) = indpost;

            good = isfinite(ind);

            %allocate the matrix
            A = NaN(size(ind));
            ind = ind(good);
            
            %and assign it the values in tevent, as ordered in ind
            tt = tevent(ch,goodevent);
            A(good) = tt(ind);
            %set it per channel
            N{1}(:,:,ch) = A;
            
            %run through the additional data, too
            for i = 2:length(N),
                tt = extra{i-1}(ch,goodevent);
                A(good) = tt(ind);
                N{i}(:,:,ch) = A;
            end;
        else
            %construct a cell matrix with all of the events in between each of the pulses
            C = mat2cell(tevent(ch,goodevent)',eventperpulse,1);
            %turn the pulse times into a cell array, too
            P = num2cell(pulseon(goodpulse))';
            
            %logical array for events that are within the requested time ranges
            if (any(rng < 0))
                ispre = cellfun(@(x,t) ((x >= t+rng(1)) & (x <= t+rng(2))), C(1:end-1),P, ...
                    'UniformOutput',false);
                pre = cellfun(@(x,c) (x(c)), C(1:end-1),ispre, ...
                    'UniformOutput',false);
            end;
            if (any(rng > 0))
                ispost = cellfun(@(x,t) ((x >= t+rng(1)) & (x <= t+rng(2))), C(2:end),P, ...
                    'UniformOutput',false);
                post = cellfun(@(x,c) (x(c)), C(2:end),ispost, ...
                    'UniformOutput',false);
            end;

            if (all(rng < 0))
                %all events are pre
                N{1}(pulseord,ch) = pre;
            elseif (all(rng >= 0))
                %all events are post;
                N{1}(pulseord,ch) = post;
            else
                %cat them
                N{1}(pulseord,ch) = cellfun(@(a,b) [a;b], pre,post, 'UniformOutput',false);
            end;

            %now run through the rest of the inputs and do the same thing
            for i = 2:length(N),
                C = mat2cell(extra{i-1}(ch,goodevent)', eventperpulse,1);

                if (any(rng < 0))
                    pre = cellfun(@(x,c) (x(c)), C(1:end-1),ispre, ...
                        'UniformOutput',false);
                end;
                if (any(rng > 0))
                    post = cellfun(@(x,c) (x(c)), C(2:end),ispost, ...
                        'UniformOutput',false);
                end;
                
                if (all(rng < 0))
                    %all events are pre
                    N{i}(pulseord,ch) = pre;
                elseif (all(rng >= 0))
                    %all events are post;
                    N{i}(pulseord,ch) = post;
                else
                    %cat them
                    N{i}(pulseord,ch) = cellfun(@(a,b) [a;b], pre,post, 'UniformOutput',false);
                end;
            end;
        end;
    end;
    
    %turn cell arrays into matrices, if necessary
    if (iscell(N{1})),
        switch opt.out,
            case 'matrixout',
                A = N{1};
                %look for non-empty elements
                good = ~cellfun(@isempty,A);
                %and get the maximum length, so that we can set up the array in advance
                len = max(cellfun(@length,A(:)));
                %allocate the array
                N{1} = NaN(len,npulse,nchan);
                %and cat non-empty elements
                N{1}(:,good) = catuneven(2,A{good});
                for i = 2:length(N),
                    A = N{i};
                    N{i} = NaN(len,npulse,nchan);
                    N{i}(:,good) = catuneven(2,A{good});
                end;
            case 'cellout',
                for i = 1:length(N),
                    N{i} = reshape(N{i},[1 npulse nchan]);
                end;
        end;
    end;
    
    varargout = N;
end;

        


