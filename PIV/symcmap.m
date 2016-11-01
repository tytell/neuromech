function varargout = symcmap(cmap0,ctr,nsteps)

if (nargin < 3),
    nsteps = Inf;
    if (nargin < 2),
        ctr = 0;
    end;
end;

if ctr > 1
    error('ctr parameter is a fraction < 1');
end

if (isnumeric(cmap0) && (size(cmap0,1) ~= 3)),
    error('cmap can only have three entries.');
end;

if (ischar(cmap0)),
    name = cmap0;
    cmap0 = zeros(3,3);
    for i = 1:length(name),
        switch lower(name(i)),
            case 'y',
                cmap0(i,:) = [1 1 0];
            case 'm',
                cmap0(i,:) = [1 0 1];
            case 'c',
                cmap0(i,:) = [0 1 1];
            case 'r',
                cmap0(i,:) = [1 0 0];
            case 'g',
                cmap0(i,:) = [0 1 0];
            case 'b',
                cmap0(i,:) = [0 0 1];
            case 'w',
                cmap0(i,:) = [1 1 1];
            case 'k',
                cmap0(i,:) = [0 0 0];
            case {'1','2','3','4','5','6','7','8','9'}
                % gray levels
                col = str2double(name(i));
                cmap0(i,:) = [col col col]/10;
        end;
    end;
end;

if (~isfinite(nsteps)),
    ctrsteps = round(ctr*256);
    nsteps2 = floor((256-ctrsteps)/2);
else
    ctrsteps = round(ctr*nsteps);
    nsteps2 = floor((nsteps+1-ctrsteps)/2);
end;

pt = linspace(0,1,nsteps2)';

cmap(1:nsteps2,:) = repmat(cmap0(2,:)-cmap0(1,:),[nsteps2 1]) .* ...
    repmat(pt,[1 3]) + repmat(cmap0(1,:),[nsteps2 1]);
cmap(end+(0:ctrsteps-1),:) = repmat(cmap0(2,:),[ctrsteps 1]);
cmap(end+(0:nsteps2-1),:) = repmat(cmap0(3,:)-cmap0(2,:),[nsteps2 1]) .* ...
    repmat(pt,[1 3]) + repmat(cmap0(2,:),[nsteps2 1]);

if (nargout == 0),
    colormap(cmap);
else
    varargout{1} = cmap;
end;



    
    
