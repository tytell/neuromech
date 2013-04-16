function [x,y,u,v,err,I] = loadInsight(files,avi,frnum)
% function [x,y,u,v,err,I] = loadInsight(files,avi,frnum)

if (nargin == 1)
	avi = [];
	frnum = [];
end;

if (ischar(files)),
	files = {files};
end;

for i = 1:length(files),
	fn = files{i};
	
	fid = fopen(fn,'r');
	if (fid ~= -1),
		title = fgets(fid);
		
        tok = regexp(title, '(\w+)=(\S+)','tokens');
        insightdata = struct();
        for i = 1:length(tok)
            num = regexp(tok{i}{2},'^"?[\d.]+"?$','once','tokens');
            if ~isempty(num)
                insightdata.(tok{i}{1}) = str2double(num);
            elseif (~isempty(tok{i}{2}) && (tok{i}{2}(1) == '"') && ...
                    (tok{i}{2}(end) == '"'))
                insightdata.(tok{i}{1}) = tok{i}{2}(2:end-1);
            else
                insightdata.(tok{i}{1}) = tok{i}{2};
            end;
        end
            
		k = findstr('I=',title);
		m = sscanf(title(k+2:end),'%f');
		k = findstr('J=',title);
		n = sscanf(title(k+2:end),'%f');
		k = findstr('Height=',title);
		h = sscanf(title(k+7:end),'%f');
		if isempty(h)
            h = 0;
        end
        
		data = fscanf(fid,'%f, %f, %f, %f, %d\n',[5 Inf]);
		fclose(fid);
		
		data = data';
		x(:,:,i) = reshape(data(:,1),[m n])';
		y(:,:,i) = h - reshape(data(:,2),[m n])';
		u(:,:,i) = reshape(data(:,3),[m n])';
		v(:,:,i) = -reshape(data(:,4),[m n])';
		err(:,:,i) = reshape(data(:,5),[m n])';
	else
		warning(sprintf('Could not open file %s. Skipping it...',fn));
	end;
end;

if (~isempty(avi)),
	if (length(frnum) ~= length(files))
		warning('Different number of images than vector files.  Ignoring images.');
	else
		mov = aviread(avi,frnum);
		
		for i = 1:length(frnum),
			I(:,:,i) = frame2im(mov(i));
		end;
	end;
end;
