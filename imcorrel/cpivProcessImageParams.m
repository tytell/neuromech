function [id,p] = cpivProcessImageParams(p)

id = [];
if (iscellstr(p{1})),
    if (length(p{1}) == 1)
        error('You must have at least two images.');
    end;

    id.FileNames = p{1};
    id.NFrames = size(id.FileNames,1)-1;

    if ((size(p{1},1) == 1) && (size(p{1},2) > 1)),
        id.FileNames = id.FileNames';
    elseif ((size(p{1},1) == 2) && (size(p{1},2) > 2)),
        id.FileNames = id.FileNames';
    end;

    if (size(id.FileNames,2) == 2),
        id.NFrames = id.NFrames + 1;
    end;

    info = imfinfo(id.FileNames{1});

    id.Width = info.Width;
    id.Height = info.Height;

    id.Frames = [1:id.NFrames; 2:id.NFrames+1]';

    i = 2;
elseif (ischar(p{1})),
    if (strfind(lower(p{1}), 'avi')),		% avi file
        id.AVIName = p{1};
        
        vid = VideoReader2(id.AVIName);
        id.NAviFrames = vid.NumberOfFrames;
        id.Width = vid.Width;
        id.Height = vid.Height;
        id.Video = vid;
        
        if (ischar(p{2})),					% second avi
            id.AVI2Name = p{2};
            
            vid2 = VideoReader2(id.AVI2Name);
            avi2len = vid2.NumberOfFrames;
            id.Video2 = vid2;
            
            if ((vid2.Width ~= id.Width) || (vid2.Height ~= id.Height)),
                error('The two AVI files must have the same size images.');
            end;
            
            i = 3;
        else
            i = 2;
        end;
        
        if (iscell(p{i})),				% defined frame skip
            a = 1;
            b = id.NAviFrames-1;
            
            frskip = p{i};
            if (length(frskip) == 1),
                dfr = frskip{1};
            elseif (length(frskip) == 2),
                a = frskip{1};
                dfr = frskip{2};
            elseif (length(frskip) == 3),
                a = frskip{1};
                dfr = frskip{2};
                b = frskip{3};
            end;
            
            id.Frames(:,1) = (a:dfr:b)';
            if (isfield(id,'AVI2Name')),
                id.Frames(:,2) = id.Frames(:,1);
            else
                id.Frames(:,2) = id.Frames(:,1)+1;
            end;
        elseif (isnumeric(p{i})),
            id.Frames = shiftdim(p{i});
            if (size(id.Frames,2) > 2),
                error('Frames must either be a vector or a two column matrix.');
            end;
            
            if (size(id.Frames,2) == 1),
                if (isfield(id,'AVI2Name')),
                    id.Frames(:,2) = id.Frames(:,1);
                else
                    id.Frames(:,2) = id.Frames(:,1)+1;
                end;
            end;

            if (any(id.Frames(:) < 1) || any(id.Frames(:) > id.NAviFrames)),
                error('Frames out of bounds.');
            end;
            if (isfield(id,'AVI2Name') && any(id.Frames(:) > avi2len)),
                error('Frames out of bounds for second AVI.');
            end;
        end;
        
        id.NFrames = size(id.Frames,1);
        i = i+1;
    else						% just an image file name
        id.I1 = im2double(imread(p{1}));
        id.Width = size(id.I1,2);
        id.Height = size(id.I1,1);
        
        id.I2 = im2double(imread(p{2}));
        if ((size(id.I2,2) ~= id.Width) || (size(id.I2,1) ~= id.Height)),
            error('Image files must be the same size');
        end;
        
        id.Frames = [1 1];
        id.NFrames = 1;
        i = 3;
    end;
elseif (isnumeric(p{1})),
    if (~isnumeric(p{2}) || (ndims(p{1}) ~= ndims(p{2})) || ...
        any(size(p{1}) ~= size(p{2}))),
        if (size(p{1},3) > 1),
            id.I1 = p{1};
            id.Frames(:,1) = (1:size(id.I1,3)-1)';
            id.Frames(:,2) = id.Frames(:,1) + 1;
            
            i = 2;
        else
            error('Images must be the same size.');
        end;
    end;

    if (~isfield(id,'I1')),
        id.I1 = p{1};
        id.I2 = p{2};
        id.Width = size(id.I1,2);
        id.Height = size(id.I1,1);
        
        id.Frames = repmat((1:size(id.I1,3))', [1 2]);
        
        i = 3;
    end;
    
    id.NFrames = size(id.Frames,1);
end;

p = p(i:end);
