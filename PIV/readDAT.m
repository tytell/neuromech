function [x,y,u,v,w, vnames] = readDAT(files,quiet)

if (nargin == 1),
    quiet = 0;
end;

if (ischar(files)),
    files = {files};
end;

if (~quiet & (length(files) > 1)),
    showwait = 1;
    timedWaitBar(0,'Loading files...');
else
    showwait = 0;
end;

for i = 1:length(files),
    fid = fopen(files{i});
    ln = fgetl(fid);                        % discard TITLE line

    ln = fgetl(fid);                        % VARIABLES line
    if (isempty(strfind(ln,'VARIABLES'))),
        fclose(fid);
        return;
    end;
    [a,b] = regexpi(ln,'"(\w+)"');
    nvar = length(a);
    for j = 1:nvar,
        vnames{j} = ln(a(j)+1:b(j)-1);
    end;

    ln = fgetl(fid);
    if (isempty(strfind(ln,'ZONE'))),
        fclose(fid);
        return;
    end;
    [a,b,tok] = regexpi(ln,'(\w+)=(\w+)');
    for j = 1:length(a),
        switch ln(tok{j}(1,1):tok{j}(1,2)),
         case 'I',
          c = str2num(ln(tok{j}(2,1):tok{j}(2,2)));
         case 'J',
          r = str2num(ln(tok{j}(2,1):tok{j}(2,2)));
        end;
    end;

    if (i == 1),
        x = zeros(r,c,length(files));
        y = zeros(r,c,length(files));
        u = zeros(r,c,length(files));
        v = zeros(r,c,length(files));
        if (nvar == 5),
            w = zeros(r,c,length(files));
        end;
    else
        if ((r ~= size(x,1)) | (c ~= size(x,2))),
            error('File %s is a different size than others.', files{i});
        end;
    end;

    if (nvar == 4),
        data = fscanf(fid,'%f %f %f %f\n',[4 Inf]);
    else
        data = fscanf(fid,'%f %f %f %f %f\n',[5 Inf]);
    end;

    fclose(fid);

    data = data';
    x(:,:,i) = reshape(data(:,1),[c r])';
    y(:,:,i) = reshape(data(:,2),[c r])';
    u(:,:,i) = reshape(data(:,3),[c r])';
    v(:,:,i) = reshape(data(:,4),[c r])';
    if (nvar == 5),
        w(:,:,i) = reshape(data(:,5),[c r])';
    end;

    if (showwait),
        if (~timedWaitBar(i/length(files))),
            break;
        end;
    end;
end;

if (showwait),
    timedWaitBar(1);
end;


