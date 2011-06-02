function camera = readCameraInfo(file)
% function camera = readCameraInfo(file)
% Reads the .cih files saved by the Photron high speed cameras
% Produces a structure called camera with various useful fields

fid = fopen(file);
done = 0;

while ~done,
    ln = fgetl(fid);

    if (ln == -1),
        done = 1;
        continue;
    end;

    if (ln(1) == '#'),
        continue;
    end;

    [a,b,tok] = regexp(ln,'^(.+) : (.+)$');
    if (isempty(a)),
        continue;
    end;

    tok = tok{1};
    thing = ln(tok(1,1):tok(1,2));
    value = ln(tok(2,1):tok(2,2));

    switch thing,
     case 'Date',
      camera.Date = value;
     case 'Time',
      camera.Time = value;
     case 'Camera Type',
      camera.Type = value;
     case 'Record Rate(fps)',
      camera.FPS = str2num(value);
     case 'Shutter Speed(s)',
      [a,b] = regexp(value,'\d+');
      camera.Shutter = str2num(value(a(1):b(1))) / ...
          str2num(value(a(2):b(2)));
     case 'Total Frame',
      camera.NFrames = str2num(value);
     case 'Start Frame',
      camera.StartFrame = str2num(value);
    end;
end;

fclose(fid);
