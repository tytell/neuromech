classdef (CaseInsensitiveProperties=true, TruncatedProperties=true) ...
         VideoReader2 < hgsetget
% VIDEOREADER2 Create a multimedia reader object that handles uncompressed AVIs.
%
%   Otherwise identical to VideoReader.


    %------------------------------------------------------------------
    % General properties (in alphabetic order)
    %------------------------------------------------------------------
    properties(GetAccess='public', SetAccess='private')
        Name            % Name of the file to be read.
        Path            % Path of the file to be read.
    end
    
    properties(GetAccess='public', SetAccess='private', Dependent)
        Duration        % Total length of file in seconds.
    end
    
    properties(GetAccess='public', SetAccess='public')
        Tag = '';       % Generic string for the user to set.
    end
    
    properties(GetAccess='public', SetAccess='private', Dependent) 
        Type            % Classname of the object.
    end
    
    properties(GetAccess='public', SetAccess='public')
        UserData        % Generic field for any user-defined data.
        isaviread
    end
    
    %------------------------------------------------------------------
    % Video properties (in alphabetic order)
    %------------------------------------------------------------------
    properties(GetAccess='public', SetAccess='private')
        BitsPerPixel    % Bits per pixel of the video data.
        FrameRate       % Frame rate of the video in frames per second.
        Height          % Height of the video frame in pixels.
        NumberOfFrames  % Total number of frames in the video stream. 
        VideoFormat     % Video format as it is represented in MATLAB.
        Width           % Width of the video frame in pixels.
    end
    
    %------------------------------------------------------------------
    % Undocumented properties
    %------------------------------------------------------------------
    properties(GetAccess='public', SetAccess='private')
        AudioCompression
        NumberOfAudioChannels
        VideoCompression
    end
    
    properties(Access='private', Hidden)
        vid
        info        
    end
    %------------------------------------------------------------------
    % Documented methods
    %------------------------------------------------------------------    
    methods(Access='public')
    
        %------------------------------------------------------------------
        % Lifetime
        %------------------------------------------------------------------
        function obj = VideoReader2(fileName, varargin)

            % If no file name provided.
            if nargin == 0
                error(message('MATLAB:audiovideo:VideoReader:noFile'));
            end

            obj.isaviread = ~isempty(which('aviread'));
            try
                obj.vid = VideoReader(fileName, varargin{:});
            catch err
                if (strcmp(err.identifier, 'MATLAB:audiovideo:VideoReader:FileInit') && ...
                        obj.isaviread)
                    obj.vid = struct([]);
                    w = warning('off','MATLAB:audiovideo:aviinfo:FunctionToBeRemoved');
                    obj.info = aviinfo(fileName);       %#ok
                    warning(w);
                    
                    fullname = VideoReader.getFullPathName(fileName);
                    [pn,fn,ext] = fileparts(fullname);
                    set(obj,'Path',pn);
                    set(obj,'Name',[fn ext]);
                else
                    rethrow(err);
                end
            end
                
            % Set properties that user passed in.
            if nargin > 1
                set(obj, varargin{:});
            end
        end


        %------------------------------------------------------------------
        % Operations
        %------------------------------------------------------------------        
        function varargout = read(obj, varargin)
            if ~isempty(obj.vid)
                v = read(obj.vid,varargin{:});
                varargout = {v};
            else
                w = warning('off','MATLAB:audiovideo:aviread:FunctionToBeRemoved');
                fr = aviread(fullfile(obj.Path,obj.Name),varargin{:});  %#ok
                warning(w);
                varargout = {fr.cdata};
            end
        end
        
        %------------------------------------------------------------------        
        % Overrides of hgsetset
        %------------------------------------------------------------------        
        function getdisp(obj)
            if ~isempty(obj.vid)
                obj.vid.getdisp();
            else
                getdisp@hgsetget(obj);
            end
        end
        function setdisp(obj)
            if ~isempty(obj.vid)
                obj.vid.setdips();
            else
                setdisp@hgsetget(obj);
            end
        end

        %------------------------------------------------------------------        
        % Overrides of builtins
        %------------------------------------------------------------------ 
        function disp(obj)
            if ~isempty(obj.vid)
                obj.vid.display();
            else
                fprintf('Uncompressed avi ''%s''\n', obj.Name);
            end
        end

        function display(obj)
            disp(obj);
        end
    end
    
    %------------------------------------------------------------------
    % Custom Getters/Setters
    %------------------------------------------------------------------
    methods
        % Properties that are not dependent on underlying object.
        function set.Tag(obj, value)
            if ~(ischar(value) || isempty(value))
                error(message('MATLAB:audiovideo:VideoReader:TagMustBeString'));
            end
            obj.Tag = value;
        end
        
        function value = get.Type(obj)
            value = class(obj);
        end
        
        % Properties that are dependent on underlying object.
        function value = get.Duration(obj)
            if ~isempty(obj.vid)
                value = obj.vid.Duration;
            else
                value = obj.info.NumFrames / obj.info.FramesPerSecond;
            end
        end
        
        function value = get.Name(obj)
            if ~isempty(obj.vid)
                value = get(obj.vid,'Name');
            else
                value = obj.Name;
            end
        end
        
        function value = get.Path(obj)
            if ~isempty(obj.vid)
                value = get(obj.vid,'Path');
            else
                value = obj.Path;
            end
        end
        
        function value = get.BitsPerPixel(obj)
            if ~isempty(obj.vid)
                value = obj.vid.BitsPerPixel;
            else
                value = log2(obj.info.NumColormapEntries);
            end
        end
        
        function value = get.FrameRate(obj)
            if ~isempty(obj.vid)
                value = obj.vid.FrameRate;
            else
                value = obj.info.FramesPerSecond;
            end
        end
        
        function value = get.Height(obj)
            if ~isempty(obj.vid)
                value = obj.vid.Height;
            else
                value = obj.info.Height;
            end
        end
        
        function value = get.NumberOfFrames(obj)
            if ~isempty(obj.vid)
                value = obj.vid.NumberOfFrames;
            else
                value = obj.info.NumFrames;
            end
        end
        
        function value = get.VideoFormat(obj)
            if ~isempty(obj.vid)
                value = obj.vid.VideoFormat;
            else
                value = 'Uncompressed AVI';
            end
        end
        
        function value = get.Width(obj)
            if ~isempty(obj.vid)
                value = obj.vid.Width;
            else
                value = obj.info.Width;
            end
        end
        
        function value = get.AudioCompression(obj)
            if ~isempty(obj.vid)
                value = obj.vid.AudioCompression;
            else
                warning('VideoReader2:NoSuchProperty',...
                    'No AudioCompression property in uncompressed AVIs');
                value = '';
            end
        end
        
        function value = get.NumberOfAudioChannels(obj)
            if ~isempty(obj.vid)
                value = obj.vid.NumberOfAudioChannels;
            else
                warning('VideoReader2:NoSuchProperty',...
                    'No NumberOfAudioChannels property in uncompressed AVIs');
                value = '';
            end
        end
        
        function value = get.VideoCompression(obj)
            if ~isempty(obj.vid)
                value = obj.vid.VideoCompression;
            else
                value = 'Uncompressed AVI';
            end
        end
    end
    
    %------------------------------------------------------------------        
    % Undocumented methods
    %------------------------------------------------------------------
    methods (Access='public', Hidden)
        
        %------------------------------------------------------------------
        % Lifetime
        %------------------------------------------------------------------
        function delete(obj)
            % Delete VideoReader object.
            if ~isempty(obj.vid)
                delete(obj.vid);
            end
        end
   
    end
end
