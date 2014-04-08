function [x,y,u,v,err,I,insightdata] = loadInsight(files,avi,frnum)
% function [x,y,u,v,err,I,insightdata] = loadInsight(files,avi,frnum)

if (nargin < 3)
    frnum = [];
    if (nargin == 1)
        avi = [];
    end;
end;

if (ischar(files)),
	files = {files};
end;

insightdata = struct([]);
for i = 1:length(files),
	fn = files{i};
	
	fid = fopen(fn,'r');
	if (fid ~= -1),
		title = fgets(fid);
		
        tok = regexp(title, '(\w+)=(([^ ,"]+)|("[^"]+"(,\s*"[^"]+")*))','tokens');
        for j = 1:length(tok)
            vals = regexp(tok{j}{2}, '"([^"]+)"', 'tokens');
            if isempty(vals)
                if ~isempty(regexp(tok{j}{2}, '^[0-9.e]+$', 'once'))
                    insightdata1.(tok{j}{1}) = str2double(tok{j}{2});
                else
                    insightdata1.(tok{j}{1}) = tok{j}{2};
                end;
            elseif ((length(vals) == 1) && ~isempty(regexp(vals{1}{1}, '^[0-9.e]+$', 'once')))
                insightdata1.(tok{j}{1}) = str2double(vals{1}{1});
            elseif (length(vals) == 1)
                insightdata1.(tok{j}{1}) = vals{1}{1};
            else
                insightdata1.(tok{j}{1}) = cellfun(@(x) x{1}, vals, 'UniformOutput',false);
            end                
        end

        if (~isfield(insightdata1,'I') || ~isfield(insightdata1,'J'))
            error('Cannot parse data file');
        end
        
        m = insightdata1.I;
        n = insightdata1.J;
        
        if (isfield(insightdata1,'Height'))
            h = insightdata1.Height;
        elseif (isfield(insightdata1,'SourceImageHeight') && strcmp(insightdata1.LengthUnit,'pixel'))
            h = insightdata1.SourceImageHeight;
        else
            h = 0;
        end
        
		data = fscanf(fid,'%f, %f, %f, %f, %d\n',[5 Inf]);
		fclose(fid);

        if (i == 1)
            x = zeros([n m length(files)]);
            y = zeros([n m length(files)]);
            u = zeros([n m length(files)]);
            v = zeros([n m length(files)]);
            err = zeros([n m length(files)]);
        end
		data = data';
		x(:,:,i) = reshape(data(:,1),[m n])';
		y(:,:,i) = h - reshape(data(:,2),[m n])';
		u(:,:,i) = reshape(data(:,3),[m n])';
		v(:,:,i) = -reshape(data(:,4),[m n])';
		err(:,:,i) = reshape(data(:,5),[m n])';
        
        insightdata = makestructarray(insightdata,insightdata1);
	else
		warning(sprintf('Could not open file %s. Skipping it...',fn));
	end;
end;

if (~isempty(avi)),
    vid = VideoReader(avi);
    if (isempty(frnum))
        frnum = 1:vid.NumberOfFrames;
    end
	if (length(frnum) ~= length(files))
		warning('Different number of images than vector files.  Ignoring images.');
    else
        I = zeros(vid.Height,vid.Width, length(frnum));
		
		for i = 1:length(frnum),
			I(:,:,i) = read(vid, frnum(i));
		end;
	end;
else
    I = [];
end

