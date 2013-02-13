function varargout = plotspikes(spiket,varargin)

if ((nargin >= 2) && ...
    isnumeric(varargin{1}) && (ndims(varargin{1}) == ndims(spiket)) && ...
        all(size(varargin{1}) == size(spiket))),
    spike = varargin{1};
    isspikeheight = true;
    i = 2;
else
    spike = ones(size(spiket));
    isspikeheight = false;
    i = 1;
end;

plottype = 'line';
spacechannels = false;
noplot = false;
lnspec = {};

while (i <= length(varargin)),
    if (matchlinespec(varargin{i})),
        lnspec = varargin(i);
        i = i+1;
    else
        switch lower(varargin{i}),
            case {'line','bars'},
                plottype = varargin{i};
                i = i+1;

            case 'spacechannels',
                spacechannels = true;
                i = i+1;

            case 'noplot',
                noplot = true;
                i = i+1;

            otherwise,
                error('Unrecognized argument %s',varargin{i});
        end;
    end;
end;

spiket = shiftdim(spiket,-1);           % make space in the first dimension
sz = ones(1,ndims(spiket));
sz(1) = 3;
spiket = repmat(spiket,sz);             % replicate the spike times three times

switch plottype,
    case 'line',
        spike = shiftdim(spike,-1);
        % y value is 0, then the spike height, then 0 again
        spike = cat(1,zeros(size(spike)),spike,zeros(size(spike)));

    case 'bars',
        spike = shiftdim(spike,-1);
        % y value is 0, then the spike height, then NaN
        spike = cat(1,zeros(size(spike)),spike,nans(size(spike)));
end;

sz = size(spiket);
sz = [sz(1)*sz(2) sz(3:end)];
if (length(sz) == 1),
    sz(2) = 1;
end;
spiket = reshape(spiket,sz);
spike = reshape(spike,sz);

if (spacechannels),
    if (ndim(spiket) > 2),
        error('Cannot space channels with more than two dimensions');
    end;
    
    if (isspikeheight),
        h = prctile(spike(:),95);
    else
        h = 1;
    end;
    
    spike = spike + h * repmat(1:size(spike,2),[size(spike,1) 1]);
end;

if (noplot),
    varargout = {spiket,spike};
else
    h = plot(spiket,spike,lnspec{:});
    
    if (nargout == 1),
        varargout = {h};
    end;
end;
