function varargout = loadPIVFilter(files, frames, fps,cutoff,ord, quiet)
% function [x,y,u,v,w?] = loadPIVFilter(files, frames, fps,cutoff,ord, quiet)

if (nargin == 5),
    quiet = false;
end;

% round the order up to the nearest odd number, if necessary
ord2 = floor(ord/2);
ord = 2*ord2 + 1;

% sort out which frames we'll actually need to load to do the filtering
frames = frames(:)';
nfiles = length(files);
nframes = length(frames);
if (ord == 1),
    filtframes = frames;
else
    filtframes = repmat(frames,[ord 1]) + repmat((-ord2:ord2)',[1 ...
                        nframes]);
end;

% eliminate out of range frame numbers
k = find((filtframes < 1) | (filtframes > nfiles));
filtframes(k) = NaN;

%check the file type
switch lower(files{1}(end-2:end)),
 case 'dat',
  filetype = 'dat';
  loadfcn = @readDAT;
  optparam = {1};

 case {'vec','vc7'},
  filetype = 'davis';
  loadfcn = @readimxvec;
  optparam = {};

 otherwise,
  error('Unrecognized file type');
end;

% load the first frame to get data sizes
[x,y,up,vp,wp,units] = feval(loadfcn,files{1},optparam{:});
nx = size(x,2);
ny = size(x,1);
if (isempty(wp)),
    isW = false;
else
    isW = true;
end;

if (ord > 1),
    %set up the filter - fir1 generates a sinc function with ord elements that
    %will cut off at a fraction of the Nyquist frequency (fps/2)
    filt = fir1(ord-1, cutoff/(fps/2));
    if (filt(1) > 0.5*filt(ord2+1)),
        warning('Filter order seems to be too low.');
    end;

    %expand it so that we can do it for each element in the frame quickly
    filt = repmat(shiftdim(filt,-1),[ny nx 1]);  
end;

% u1 etc are 3D matrices with the 3rd dimension being the number of elements
% we need for filtering
prevframes = zeros(ord,1);
u1 = zeros(ny,nx,ord);
v1 = zeros(ny,nx,ord);
valid = true(ny,nx,ord);
u = repmat(NaN,[ny nx nframes]);
v = repmat(NaN,[ny nx nframes]);
if (isW),
    w1 = zeros(ny,nx,ord);
    w = repmat(NaN,[ny nx nframes]);
end;

if (~quiet),
    timedWaitBar(0,'Filtering...');
end;

for i = 1:length(frames),
    if (all(isnan(filtframes(:,i)))),
        warning('No frames available to get frame %d.',i);
        continue;
    end;

    % check which frames we loaded before
    [c,prev,common] = intersect(prevframes,filtframes(:,i));

    % we'll need them in a new position, though (i.e., we had frames
    % [1 2 3 4] and now we want [3 4 5 6], so we move 3 and 4 to the
    % new position, indicated by common
    u1(:,:,common) = u1(:,:,prev);
    v1(:,:,common) = v1(:,:,prev);

    % fill frames that we don't have with zeros
    bad = find(isnan(filtframes(:,i)));
    u1(:,:,bad) = 0;
    v1(:,:,bad) = 0;

    if (isW),
        w1(:,:,common) = w1(:,:,prev);
        w1(:,:,bad) = 0;
    end;

    valid(:,:,common) = valid(:,:,prev);
    valid(:,:,bad) = false;
    
    % and remove the bad frames from the new ones we'll need to load
    new = setdiff(1:ord,common);
    new = setdiff(new,bad);

    % load the new frames
    if (~isempty(new)),
        for j = 1:length(new),
            if (isW),
                [x1,y1,u1(:,:,new(j)),v1(:,:,new(j)),w1(:,:,new(j))] = ...
                    feval(loadfcn,files{filtframes(new(j),i)},optparam{:});
            else
                [x1,y1,u1(:,:,new(j)),v1(:,:,new(j))] = ...
                    feval(loadfcn,files{filtframes(new(j),i)},optparam{:});
            end;
        end;
    else
        fprintf('Step %d: No new frames.\n',i);
    end;

    if (strcmp(filetype,'dat')),
        %DAT files have zeros for missing elements, so make sure we
        %only sum over the valid (non-zero) elements
        valid(:,:,new) = (u1(:,:,new) ~= 0) | (v1(:,:,new) ~= 0);
        if (isW),
            valid(:,:,new) = valid(:,:,new) | (w1(:,:,new) ~= 0);
        end;
    else
        valid(:,:,new) = isfinite(u1(:,:,new)) & isfinite(v1(:,:,new));
        if (isW),
            valid(:,:,new) = valid(:,:,new) & isfinite(w1(:,:,new));
        end;
    end;
    
    if (ord > 1),
        u1(~valid) = 0;
        v1(~valid) = 0;
        if (isW),
            w1(~valid) = 0;
        end;
        valid = double(valid);
        
        %We'll need to divide by the sum of the filter, which is 1 if we
        %have all the elements, but since some may be invalid, we have to
        %sum it each time
        denom = sum(filt.*valid, 3);
        denom(denom == 0) = 1;

        %Do the filter
        u2 = sum(u1.*filt,3)./denom;
        v2 = sum(v1.*filt,3)./denom;
        %return elements that were NaN in the current frame back to NaN
        u2(~valid(:,:,ord2+1)) = NaN;
        v2(~valid(:,:,ord2+1)) = NaN;
        u(:,:,i) = u2;
        v(:,:,i) = v2;
        if (isW),
            w2 = sum(w1.*filt,3)./denom;
            w2(~valid(:,:,ord2+1)) = NaN;
            w(:,:,i) = w2;
        end;
    else
        u(:,:,i) = u1;
        v(:,:,i) = v1;
        if (isW),
            w(:,:,i) = w1;
        end;
    end;
           
    prevframes = filtframes(:,i);

    if (~quiet),
        timedWaitBar(i/nframes);
    end;
end;

if (~quiet),
    timedWaitBar(1);
end;

if (isW),
    varargout = {x,y,u,v,w};
else
    varargout = {x,y,u,v};
    if (nargout == 5),
        varargout{5} = [];
    end;
end;

    
    