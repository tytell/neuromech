function varargout = identifyVortices(x,y,u,v, varargin)
% function vxr = identifyVortices(x,y,u,v, ...)
%                identifyVortices(...,'thresh',thresholdDCEV)
%                identifyVortices(...,'minframes',minimumFrameDuration)
%                identifyVortices(...,'mincirc',minimumCirculation)
%                identifyVortices(...,'quiet')
%                identifyVortices(...,'circmagnify',radiusIncrease)
%                identifyVortices(...,'nocirc')
%     or   [vxr,label] = indentify(...)
%
% thresh is the DCEV threshold for finding vortices.  Err on the side
% of too low.
% minframes is the minimum number of frames we need to find a vortex
% in in order to keep it
% mincirc is the minimum circulation magnitude a vortex must have, in 
% at least one frame in which it's defined, to be kept.  NB: can't have
% nocirc option if mincirc > 0
% circmag defines how much we increase the geometric radius to
% estimate circulation

thresh = 300;
minframes = 1;
quiet = false;
mincirc = 0;
%increase the radius slightly for circulation calculations, because we
%step out from the center and find the max
circmag = 1.5;
calccirc = true;

opts = varargin;
i = 1;
while (i <= length(opts)),
    switch lower(opts{i}),
     case {'thresh','threshold'},
      thresh = opts{i+1};
      i = i+2;
     case 'minframes',
      minframes = opts{i+1};
      i = i+2;
     case 'quiet',
      quiet = true;
      i = i+1;
     case 'circmagnify',
      circmag = opts{i+1};
      i = i+2;
     case 'mincirc',
      mincirc = opts{i+1};
      i = i+2;
     case 'nocirc',
      calccirc = false;
      i = i+1;
     case 'dcev',
     %we already had dcev, so don't recalculate it here
      dcev = opts{i+1};
      i = i+2;
     otherwise,
      warning('Unrecognized option %s.',opts{i});
      i = i+1;
    end;
end;

if ((mincirc > 0) & ~calccirc),
    warning(['Circulation must be calculated if a minimum circulation is ' ...
             'specified.']);
    calccirc = true;
end;

if (nargout == 2),
    returnlab = true;
else
    returnlab = false;
end;

if (~quiet),
    timedWaitBar(0,'Calculating DCEV...');
end;

%construct a "signed dcev" (dcev signed by the vorticity at that point)
if (~exist('dcev')),
    DU = velderiv(x,y,u,v);
    dcev = disccomplexeig(DU);
    vort = cat(3,DU.dvdx) - cat(3,DU.dudy);
    dcev = dcev .* sign(vort);
end;

if (~quiet),
    timedWaitBar(0,'Thresholding...');
end;

%threshold the data
Dplus = dcev > thresh;
Dminus = dcev < -thresh;

%return something reasonable if we don't find any vortices
%usually this comes from a units problem (e.g., position in mm,
%and velocity in m/s)
if (all(Dplus(:) == 0) & all(Dminus(:) == 0)),
    warning('No vortices found.  Are units right?');

    vxr = struct('ctrx',{},'ctry',{},'rgeom',{},...
                 'majordiam',{},'minordiam',{},'angle',{},...
                 'rcirc',{},'circmax',{},'frames',{});
    varargout{1} = vxr;

    if (returnlab),
        varargout{2} = zeros(size(Dplus));
    end;
    return;
end;

%18-connected neighborhood for semi-conservative detection
if (~quiet),
    timedWaitBar(0,'Finding clockwise regions...');
end;
L = bwlabeln(Dplus,18);
if (~quiet),
    timedWaitBar(0,'Finding counterclockwise regions...');
end;
Lminus = bwlabeln(Dminus,18);
Nplus = max(L(:));
L = L + (Lminus ~= 0).*(Lminus + Nplus);

if (~quiet),
    timedWaitBar(0,'Checking region duration...');
end;

%eliminate regions that span too few frames
%first get the frame ranges
startframe = (size(L,3)+1)*ones(1,max(L(:)));
endframe = zeros(1,max(L(:)));
for i = 1:size(L,3),
    ind = L(:,:,i);         % find the regions present in this frame
    ind = ind(ind ~= 0);    % ind is a list of regions in this frame

    %NB: this is a little tricky.  ind will contain the index of the
    %vortex, but repeated many times (for each pixel in an identified 
    %vortex).
    %but there's no need to use unique or something like that - we
    %can just use the index, with repeats and all, in the assign.
    startframe(ind) = min(cat(1,i*ones(1,length(ind)), ...
                              startframe(ind)));
    endframe(ind) = max(cat(1,i*ones(1,length(ind)), ...
                              endframe(ind)));
end;

if (~quiet),
    timedWaitBar(0,'Sorting vortices...');
end;

dur = endframe - startframe + 1;

%sort vortices by starting frame (ascending) and duration (descending)
good = find(dur >= minframes);
[q,ord] = sortrows([startframe(good)' (-dur(good)')]);
ord = good(ord);

if (~quiet),
    fprintf('%d of %d (%d%%) eliminated because of short durations.\n',...
            length(dur)-length(good), length(dur), ...
            round((length(dur)-length(good))/length(dur)*100));
end;

%and renumber them to reflect the ones that are too short
L2 = zeros(size(L));
szxy = size(L2,1)*size(L2,2);

%set up the vxr structure
empty = cell(length(good),1);
vxr = struct('ctrx',empty,'ctry',empty,'rgeom',empty,...
             'majordiam',empty,'minordiam',empty,'angle',empty,...
             'rcirc',empty,'circmax',empty,'frames',empty);

for i = 1:length(ord),
    %get the frame numbers
    fr1 = startframe(ord(i)):endframe(ord(i));
    vxr(i).frames = fr1;

    %find the vortex number in the original label matrix
    ind = find(L(:,:,fr1) == ord(i));
    ind = ind + (fr1(1)-1)*szxy;

    %and replace it in the new label matrix with the right
    %order number
    L2(ind) = i;
end;
L = L2;

theta = linspace(0,2*pi,33);

indx = 1:size(x,2);
indy = 1:size(y,1);

if (~quiet),
    timedWaitBar(0,'Estimating boundaries...');
end;

%identify circles for each region and calculate circulation
for i = 1:size(u,3),
    %only look at those regions that show up on this frame
    L1 = L(:,:,i);
    rgns = unique(L1(:));
    rgns = rgns(2:end);

    %if we didn't find any regions, then just skip to the next frame
    if (isempty(rgns)),
        continue;
    end;

    L1 = L1 - rgns(1) + 1;
    L1(L1 < 0) = 0;
    %we subtract off the lowest region number in this frame, because
    %regionprops creates a structure array with the true numbers of the
    %regions.  So if the lowest region number in this frame is 10000, but
    %there are only 4 regions, regionprops would normally create a length
    %10004 structure
    
    %collect the appropriate parameters
    stats = regionprops(L1,'Centroid','MajorAxisLength','MinorAxisLength',...
                        'Orientation');

    for j = 1:length(rgns),
        a = rgns(j);
        b = rgns(j) - rgns(1) + 1;

        %Centroid is a (potentially) fractional index value - interpolate
        %the correct x and y values based on the x and y matrices
        vxr(a).ctrx(end+1) = interp1(indx,x(1,:),stats(b).Centroid(1));
        vxr(a).ctry(end+1) = interp1(indy,y(:,1),stats(b).Centroid(2));

        %rgeom is the average geometric radius, based on the thresholding
        vxr(a).rgeom(end+1) = sqrt(stats(b).MajorAxisLength * ...
                                   stats(b).MinorAxisLength / 4);
        vxr(a).majordiam(end+1) = stats(b).MajorAxisLength;
        vxr(a).minordiam(end+1) = stats(b).MinorAxisLength;
        vxr(a).angle(end+1) = stats(b).Orientation*pi/180;

        if (calccirc),
            %calculate circulation, stepping out from the middle
            [circ,circdetails] = circulation(x,y,u(:,:,i),v(:,:,i), ...
                                       vxr(a).ctrx(end),vxr(a).ctry(end), ...
                                       [circmag*vxr(a).majordiam(end)/2 ...
                                        circmag*vxr(a).minordiam(end)/2], ...
                                       vxr(a).angle(end));
            %record the max value
            vxr(a).circmax(end+1) = circ;

            %and the radius where we found it
            vxr(a).rcirc(end+1) = circdetails.rAtMaxCirc;
        end;
    end;

    if (~quiet),
        timedWaitBar(i/(size(u,3)+2));
    end;
end;

if (~quiet),
    timedWaitBar(1);
end;

%eliminate low circulation vortices
if (mincirc > 0),
    maxcirc = zeros(size(vxr));
    for i = 1:length(vxr),
        maxcirc(i) = max(abs(vxr(i).circmax));
    end;

    good = find(maxcirc >= mincirc);

    if (~quiet),
        fprintf('%d of %d (%d%%) eliminated because of low circulation.\n',...
                length(vxr)-length(good), length(vxr), ...
                round((length(vxr)-length(good))/length(vxr)*100));
    end;

    vxr = vxr(good);

    %if we need to return the label matrix, then we need to renumber
    %it to reflect the vortices we just eliminated
    if (returnlab),
        for i = 1:length(good),
            if (i ~= good(i)),
                %get the frame numbers
                fr1 = vxr(i).frames;

                %find the original vortex number in the label matrix
                ind = find(L(:,:,fr1) == good(i));
                ind = ind + (fr1(1)-1)*szxy;

                %and replace it in the new label matrix with the right
                %order number
                L(ind) = i;
            end;
        end;
    end;
end;

%if there are two outputs, process the label matrix so it reflects the
%number and order of vortices as we just re-ordered them
if (returnlab),
    varargout{2} = L;
end;

varargout{1} = vxr;

    




    