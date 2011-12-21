function [I1,I2] = cpivGetImages(id, fr)

fr1 = id.Frames(fr,1);
fr2 = id.Frames(fr,2);

if (isfield(id,'FileNames')),
    I1 = im2double(imread(id.FileNames{fr1,1}));
    if (size(id.FileNames,2) == 2),
        I2 = im2double(imread(id.FileNames{fr1,2}));
    else
        I2 = im2double(imread(id.FileNames{fr2,1}));
    end;
elseif (isfield(id,'Video')),
    if (isfield(id,'Video2')),
        I1 = read(id.Video,fr1);
        I2 = read(id.Video2,fr2);
    else
        I1 = read(id.Video, fr1);
        I2 = read(id.Video, fr2);
    end;
    if (size(I1,3) == 3)
        I1 = rgb2gray(I1);
        I2 = rgb2gray(I2);
    end;
elseif (isfield(id,'AVIName')),
    if (isfield(id,'AVI2Name')),
        I1 = im2double(frame2im(aviread(id.AVIName, fr1)));
        I2 = im2double(frame2im(aviread(id.AVI2Name, fr2)));
    else
        mov = aviread(id.AVIName, [fr1 fr2]);
        I1 = im2double(frame2im(mov(1)));
        I2 = im2double(frame2im(mov(2)));
    end;
elseif (isfield(id,'I1')),
    I1 = id.I1(:,:,fr1);
    I2 = id.I2(:,:,fr2);
end;

		