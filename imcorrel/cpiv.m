function varargout = cpiv(varargin)
% [x,y,u,v,data] = cpiv(Images,[OutputFile],Grid,...)
%
% The images can be specified in a variety of ways:
%     Two matrices -          cpiv(I1,I2,...)
%         Note that these matrices can be three dimensional to specify
% 	many images.
%     Two image files -       cpiv(File1,File2,...)
%     Many image files -      cpiv(FileNames,...)
%       FileNames is a cell string of file names.  Normally it's a
%       vector, in which case the ith file is correlated with the i+1th
%       file.  It can be a two column matrix, though, in which case the
%       ith row specifies the image pair for the ith vector field.
%     An AVI file -           cpiv(AVIFile,frames,...)
%       The frames specification can either be a vector, in which case
%       each frame listed in in frames is compared to the frame following
%       it in the avi (i.e., if frames = [1 6], then frame 1 will be
%       compared to frame 2 and frame 6 to frame 7), or a two column
%       matrix, where each row specifies the two frames that will be
%       compared.
%     Two AVI files -         cpiv(AVI1,AVI2,frames,...)
%       If frames is a vector, the same frame will be compared in each
%       AVI (i.e., if frames = [1 6], then frame 1 in AVI1 will be
%       compared to frame 1 in AVI2, etc.).  Frames can also be a two-
%       column matrix, as above.
%
%   The OutputFile is optional, but if specified it will be a Matlab data file
%   containing x,y,u,v, the data structure (described below) and all the
%   image and grid specifications.
%
%   The grid is specified in two ways:
%     grid = [win search offset]
%       Square interrogation regions.
%     grid = [winx winy searchx searchy offsetx offsety]
%       Rectangular interrogation regions.
%     TODO: add arbitrary specification, like in superpiv.
%
%     Any of these specifications can have multiple rows, which specifies
%     a multiple pass PIV.
%
%   You can also pass a variety of options:
%     'InterpolationMode', ['agw' or 'sts']
%       Uses adaptive gaussian window (agw) or a thin-plate smoothing
%       spline (sts) to interpolate vectors between passes.  AGW is faster
%       but less accurate and STS is slower (much slower if there are
%       many vectors), but more accurate.  See Spedding and Rignot 98.
%       Uses default agw if the number of vectors > 1000, otherwise uses
%       sts.
%     'HartCorrection', level (default = 1.0)
%       Uses Doug Hart's correction method for regions with signal to
%       noise ratio < level.
%     'SaveMode', level (default = 1)
%       Saves more data with increasing levels.  Each level saves
%       everything at the previous levels.  All of the information from 
%       levels > 0 is saved in a structure called data.
%         0 - Save only x,y,u, and v
%         1 - Save final pass signal to noise ratio (SNR), signal 
%           intensity (Signal), noise intensity (Noise), whether the 
%           Hart correction was used (isCorrected), and error values
%           (Error).
%         2 - Save intermediate pass velocities, signal, and error data.
%           All the elements have the same name as in level 1, but
%           are now cell arrays with as many elements as passes.
%         3 - Save final correlation matrix (Correlation) at each 
%           interrogation region (thus Correlation is a 4D matrix).
%         4 - Save intermediate correlation matrices.  Correlation is
%           then a cell array, also.
%     'BackgroundVelocity', [u v] (default [0 0])
%       Assumes a background flow (as in a flow tank) of (u,v) and shifts
%       the second image of each frame accordingly.  If the background flow
%       is large, this allows much smaller search regions to be used than
%       is otherwise possible.  See Spedding Exp.Fluids 03 for a similar idea.
%       Note that integral background flows are *much* faster to process
%       than fractional, and produce more or less the same results.
%     'Region', [x1 x2 y1 y2] (default whole image)
%       Only processes the passed region, specified as in axis.
%
% Returns
%   x,y - Pixel coordinates of the centers of the PIV interrogation regions.
%     NB: These are evenly spaced, but the actual PIV vectors should not be.
%     The best estimate of the position of a vector (u,v) estimated in a
%     region centered at (x,y) is actually (x+u/2,y+v/2).  cpiv takes this
%     into account when interpolating for multiple pass PIV, but not in the
%     final returned values.
%   u,v - Pixel displacements.  In a multi-frame data set, these will be 3D,
%     with the number of frames in the third dimension.
%   data - Structure containing the following elements:
%     Signal - Correlation signal intensity (from 0 to 1)
%     Noise - Height of the next highest non-signal peak in the correlation
%       matrix (also 0 to 1)
%     SNR - Signal to noise ratio.  Ratio of the above.  For some reason it
%       tends to be higher from cpiv than from other PIV routines.
%     Error - Error value.  0 is no error.  Error codes are created by
%       adding together the different error values:
%         1 - Correlation peak off the left side
%         2 - Correlation peak off the right side
%         4 - Correlation peak off the top edge
%         8 - Correlation peak off the bottom edge
%         16 - No correlation (usually means the image was all zeros in
%           that region)
%         32 - Hart correction error.  Could not apply the Hart correction
%       You can call cpivParseError with an error matrix and it will parse
%       all the different combinations of these codes and display what
%       they mean.
%     IsCorrected - 1 if the Hart correction was applied, 0 if it wasn't
%     Correlation - The actual correlation matrix, with size [r c m n], where
%       there are m x n interrogation regions, each with an r x c correlation
%       matrix.  Only included if SaveMode is 2 or more.
%
%   Additionally, when an output file is specified, cpiv creates the following
%   variables:
%    
%     imageData - Structure specifying the images or AVI files
%     gridData - Structure specifying the size of the grid
%     options - All the options passed to cpiv
%
% See Also
%   cpivParseError, quiverc, multiplot
%
% Requires
%        cpiv*.m, plus agw.m

% Deal with canceling from the wait dialog
global isCanceled;

if ((nargin == 1) & ischar(varargin{1}) & strcmp(varargin{1},'cancel')),
    isCanceled = 1;
    delete(gcf);
    return;
end;

% Don't cancel if we're just starting
isCanceled = 0;

% Process the parameters
[imageData,outputFile,gridData, options] = cpivParams(varargin);
if (~isstruct(options) && isempty(options))
    vargout = cell(1,nargout);
    return;
end;

if (isempty(outputFile) & (nargout == 0)),
    error('No output specified.');
end;

% Only track processing time if we're doing more than 1 frame
if (imageData.NFrames > 1),
    start = cputime;
    passtime = zeros(1,gridData.NPasses);
    
    hBar = waitbar(0,'Calculating...','CreateCancelBtn','cpiv(''cancel'')');
end;

for fr = 1:imageData.NFrames,
    % Get the current images from this frame
    [I1,I2] = cpivGetImages(imageData, fr);
    
    % process background flow if it exists
    if (~isempty(options.I2tform)),
        if (isstruct(options.I2tform)),
            % fractional flow requires interpolation, so use a sophisticated,
            % but slow method
            I2 = imtransform(I2,options.I2tform);
        else
            % integral displacement is easy -- we can just shift indices
            i = max(options.I2tform(2)+1,1):min(options.I2tform(2)+size(I2,1),size(I2,1));
            j = max(options.I2tform(1)+1,1):min(options.I2tform(1)+size(I2,2),size(I2,2));
            I2(i,j) = I2(i-options.I2tform(2), j-options.I2tform(1));
            
            % put zeros on the edges where necessary
            I2(1:options.I2tform(2),:) = 0;
            I2(end+options.I2tform(2)+1:end,:) = 0;
            I2(:,1:options.I2tform(1)) = 0;
            I2(:,end+options.I2tform(1)+1:end) = 0;
        end;
    end;
    
    clear data1;
    for i = 1:gridData.NPasses,
        passStart = cputime;		% track how long each pass takes
        if (fr > 1),				% but only update the waitbar after we've processed one frame
            remain = (imageData.NFrames - fr)*sum(passtime) + sum(passtime(i:end));
            done = cputime - start;
            waitbar(done/(done+remain), hBar, sprintf('Calculating frame %d (%d:%02d remaining)...', fr,floor(remain/60),ceil(mod(remain,60))));
            if (isCanceled),
                return;
            end;
        end;

        % calculate the points based on our grid specification
        [x0,y0, wsz,ssz] = cpivGetPoints(gridData, fr, i);
        x = x0;
        y = y0;
        
        % if we have a previous pass, interpolate velocities
        if (i > 1),
            [ug,vg] = cpivInterpolateVelocities(prevx,prevy,u1,v1,disperr, x,y, ...
                                                gridData,fr, options.InterpolationMode);
            
            x(:,:,2) = x + round(ug);
            y(:,:,2) = y + round(vg);
            
            % check to make sure none of our positions are off the image
            [x,y] = cpivCheckGrid(imageData.Width,imageData.Height, ssz(1),ssz(2), x,y);
        else
            ug = [];
            vg = [];
        end;
        
        % get the correlation matrix
        [phi, correrr] = imcorrel(I1,I2, x,y, wsz,ssz);

        if (options.Symmetric),
            if (size(x,3) == 2),
                x = x(:,:,[2 1]);
                y = y(:,:,[2 1]);
            end;
            
            [phi2, correrr2] = imcorrel(I2,I1, x,y, wsz,ssz);
            
            phi2 = phi2(end:-1:1,end:-1:1,:,:);
            
            if (options.SaveMode >= 3),
                phi0 = phi;
            end;
            phi = phi2.*phi;
            correrr = bitand(correrr, correrr2);
        end;
        
        % calculate displacement and signal to noise ratios
        [idx,idy, u1,v1, signal,disperr] = cpivGetDisplacement(phi);
        [snr,noise] = cpivGetSNR(phi, idx,idy, signal);
        
        % correct the correlation matrix if ratios are low
        [phiCorr, iscorrind, harterr] = cpivCorrectPhi(phi, u1,v1, snr, options.HartCorrection(i));
        iscorr = zeros(size(u1));
        iscorr(iscorrind) = 1;
        disperr = bitor(disperr, harterr*32);
        
        % and recalculate displacement and SNR for corrected values
        if (~isempty(iscorrind)),
            [idxCorr,idyCorr, uCorr,vCorr, signalCorr,disperrCorr] = cpivGetDisplacement(phiCorr, iscorrind);
            [snrCorr,noiseCorr] = cpivGetSNR(phiCorr, idxCorr,idyCorr, signalCorr, iscorrind);
            
            % insert the corrected values
            u1(iscorrind) = uCorr;
            v1(iscorrind) = vCorr;
            signal(iscorrind) = signalCorr;
            noise(iscorrind) = noiseCorr;
            disperr(iscorrind) = disperrCorr;
            snr(iscorrind) = snrCorr;
        end;
        
        % NB: if we had a previous pass, then displacements calculated 
        % are offsets from the *displaced* positions based on our guess
        % for velocity.  So we have to add the guesses back in to get a
        % real value
        if (~isempty(ug)),
            u1 = u1 + round(ug);
            v1 = v1 + round(vg);
        end;
        
        if (options.SaveMode >= 2),		% save interpolation values
            data1.x{i} = x0;
            data1.y{i} = y0;
            data1.U{i} = u1 + options.UVBkgnd(1);
            data1.V{i} = v1 + options.UVBkgnd(2);
            data1.Ug{i} = ug;
            data1.Vg{i} = vg;
            
            data1.SNR{i} = snr;
            data1.Signal{i} = signal;
            data1.Noise{i} = noise;
            data1.Error{i} = disperr;
            data1.IsCorrected{i} = iscorr;
        end;
        if (options.SaveMode >= 4),		% save all phi values
            data1.Correlation{i} = phiCorr;
            if (options.Symmetric),
                data1.CorrFwd{i} = phi0;
                data1.CorrRev{i} = phi2;
            end;
        end;
        
        % save the previous positions
        prevx = x(:,:,1);
        prevy = y(:,:,1);
        
        if (imageData.NFrames > 1),
            passtime(i) = (passtime(i)*(fr-1) + cputime - passStart)/fr;
        end;
    end;
    
    if (options.SaveMode == 1),
        data1.SNR = snr;
        data1.Signal = signal;
        data1.Noise = noise;
        data1.Error = disperr;
        data1.IsCorrected = iscorr;
    end;
    
    if (options.SaveMode == 3),		% save last phi value
        data1.Correlation = phiCorr;
        if (options.Symmetric),
            data1.CorrFwd = phi0;
            data1.CorrRev = phi2;
        end;
    end;

    if (gridData.NPasses == 1),
        cc = struct2cell(data1);
        fn = fieldnames(data1);
        for q = 1:length(cc),
            if ((length(cc{q}) == 1) & iscell(cc{q})),
                cc{q} = cc{q}{:};
            end;
        end;

        data1 = cell2struct(cc,fn,1);
    end;

    u(:,:,fr) = u1 + options.UVBkgnd(1);
    v(:,:,fr) = v1 + options.UVBkgnd(2);
    data(fr) = data1;
    
    x = x0;
    y = y0;
    
    if (~isempty(outputFile)),
        save(outputFile,'x','y','u','v', 'imageData','gridData','options');
        if (options.SaveMode > 0),
            save(outputFile,'data','-append');
        end;
    end;
end;

if (imageData.NFrames > 1),
    delete(hBar);
end;

if (nargout >= 4),
    varargout = {x(:,:,1) y(:,:,1) u v};
end;
if ((nargout == 5) & (options.SaveMode > 0)),
    varargout{5} = data;
end;

