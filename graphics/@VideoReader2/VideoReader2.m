classdef (CaseInsensitiveProperties=true, TruncatedProperties=true) ...
         VideoReader2 < VideoReader
% VIDEOREADER2 Create a multimedia reader object that handles uncompressed AVIs.
%
%   Otherwise identical to VideoReader.

%{
    %------------------------------------------------------------------
    % General properties (in alphabetic order)
    %------------------------------------------------------------------
    properties(GetAccess='public', SetAccess='private', Dependent)
        Duration        % Total length of file in seconds.
        Name            % Name of the file to be read.
        Path            % Path of the file to be read.
    end
    
    properties(GetAccess='public', SetAccess='public')
        Tag = '';       % Generic string for the user to set.
    end
    
    properties(GetAccess='public', SetAccess='private', Dependent) 
        Type            % Classname of the object.
    end
    
    properties(GetAccess='public', SetAccess='public')
        UserData        % Generic field for any user-defined data.
    end
    
    %------------------------------------------------------------------
    % Video properties (in alphabetic order)
    %------------------------------------------------------------------
    properties(GetAccess='public', SetAccess='private', Dependent)
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
    properties(GetAccess='public', SetAccess='private', Dependent, Hidden)
        AudioCompression
        NumberOfAudioChannels
        VideoCompression
    end
    
    %------------------------------------------------------------------
    % Private properties
    %------------------------------------------------------------------
    properties(Access='private', Hidden)
        % To help support future forward compatibility.
        SchemaVersion = 7.11;
        
        % To handle construction on load.
        ConstructorArgs
        
        % Enable frame counting
        % The default value of this property is TRUE.
        % If the value is FALSE, then the number of frames in the video
        % file is not determined and reported. This enables faster object
        % construction. The value held by the NumberOfFrames property of
        % the object is not valid.
        EnableFrameCounting;
    end
    
    properties(Access='private', Hidden, Transient)
        % Underlying implementation object.
        VideoReaderImpl 
    end
    %}
    
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

            try
                obj = obj@VideoReader(fileName, varargin{:});
            catch err
                
        end
    end

%{
        %------------------------------------------------------------------
        % Operations
        %------------------------------------------------------------------        
        varargout = read(obj, varargin)
        inspect(obj)
        
        %------------------------------------------------------------------        
        % Overrides of hgsetset
        %------------------------------------------------------------------        
        getdisp(obj)
        setdisp(obj)

        %------------------------------------------------------------------        
        % Overrides of builtins
        %------------------------------------------------------------------ 
        disp(obj)
        display(obj)
        c = horzcat(varargin)
        c = vertcat(varargin)
    end
    
    methods(Static)
        
        %------------------------------------------------------------------
        % Operations
        %------------------------------------------------------------------
        
        function formats = getFileFormats()
            % GETFILEFORMATS
            %
            %    FORMATS = VIDEOREADER.GETFILEFORMATS() returns an object array of 
            %    audiovideo.FileFormatInfo objects which are the formats 
            %    VIDEOREADER is known to support on the current platform. 
            %
            %    The properties of an audiovideo.FileFormatInfo object are:
            %
            %    Extension   - The file extension for this file format
            %    Description - A text description of the file format
            %    ContainsVideo - The File Format can hold video data
            %    ContainsAudio - The File Format can hold audio data
            %
            extensions = audiovideo.mmreader.getSupportedFormats();
            formats = audiovideo.FileFormatInfo.empty();
            for ii=1:length(extensions)
                formats(ii) = audiovideo.FileFormatInfo( extensions{ii}, ...
                                                         VideoReader.translateDescToLocale(extensions{ii}), ...
                                                         true, ...
                                                         false );
            end
            
            
            % sort file extension
            [~, sortedIndex] = sort({formats.Extension});
            formats = formats(sortedIndex);
            
            
        end
    end

    methods(Static, Hidden)
        %------------------------------------------------------------------
        % Persistence
        %------------------------------------------------------------------        
        obj = loadobj(B)
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
        function set.Type(obj, value)
            obj.setImplValue('Type', value);
        end
        
        % Properties that are dependent on underlying object.
        function value = get.Duration(obj)
            value = obj.getImplValue('Duration');
        end
        function set.Duration(obj, value)
            obj.setImplValue('Duration', value);
        end
        
        function value = get.Name(obj)
            value = obj.getImplValue('Name');
        end
        function set.Name(obj, value)
            obj.setImplValue('Name', value);
        end
        
        function value = get.Path(obj)
            value = obj.getImplValue('Path');
        end
        function set.Path(obj, value)
            obj.setImplValue('Path', value);
        end
        
        function value = get.BitsPerPixel(obj)
            value = obj.getImplValue('BitsPerPixel');
        end
        function set.BitsPerPixel(obj, value)
            obj.setImplValue('BitsPerPixel', value);
        end
        
        function value = get.FrameRate(obj)
            value = obj.getImplValue('FrameRate');
        end
        function set.FrameRate(obj, value)
            obj.setImplValue('FrameRate', value);
        end
        
        function value = get.Height(obj)
            value = obj.getImplValue('Height');
        end
        function set.Height(obj, value)
            obj.setImplValue('Height', value);
        end
        
        function value = get.NumberOfFrames(obj)
            value = obj.getImplValue('NumberOfFrames');
            % value = obj.NumberOfFrames;
        end
        function set.NumberOfFrames(obj, value)
            obj.setImplValue('NumberOfFrames', value);
            % obj.NumberOfFrames = value;
        end
        
        function value = get.VideoFormat(obj)
            value = obj.getImplValue('VideoFormat');
        end
        
        function set.VideoFormat(obj, value)
            obj.setImplValue('VideoFormat', value);
        end
        
        function value = get.Width(obj)
            value = obj.getImplValue('Width');
        end
        function set.Width(obj, value)
            obj.setImplValue('Width', value);
        end
        
        function value = get.AudioCompression(obj)
            value = obj.getImplValue('AudioCompression');
        end
        function set.AudioCompression(obj, value)
            obj.setImplValue('AudioCompression', value);
        end
        
        function value = get.NumberOfAudioChannels(obj)
            value = obj.getImplValue('NumberOfAudioChannels');
        end
        function set.NumberOfAudioChannels(obj, value)
            obj.setImplValue('NumberOfAudioChannels', value);
        end
        
        function value = get.VideoCompression(obj)
            value = obj.getImplValue('VideoCompression');
        end
        function set.VideoCompression(obj, value)
            obj.setImplValue('VideoCompression', value);
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
            try
                delete(obj.getImpl());
            catch exception
                VideoReader.handleImplException( exception );
            end
        end
   
        %------------------------------------------------------------------
        % Operations
        %------------------------------------------------------------------
        function result = hasAudio(obj)
            try
                result = hasAudio(obj.getImpl());
            catch exception
                VideoReader.handleImplException( exception );
            end
        end
        
        function result = hasVideo(obj)
            try
                result = hasVideo(obj.getImpl());
            catch exception 
                VideoReader.handleImplException( exception );
            end
        end
        
        function populateNumFrames(obj)
            try
                populateNumFrames(obj.getImpl());
            catch exception 
                VideoReader.handleImplException( exception );
            end
        end
    end
    
    methods (Static, Access='public', Hidden)
        
        function handleImplException(implException)  
            messageArgs = { implException.identifier };
            if (~isempty(implException.message))
                messageArgs{end+1} = implException.message;
            end
            
            msgObj = message(messageArgs{:});
            throwAsCaller(MException(implException.identifier, msgObj.getString)); 
        end
        
    end
    
    methods (Static, Access='private', Hidden)
        function errorIfImageFormat( fileName )
            isImageFormat = false;
            try 
                % see if imfinfo recognizes this file as an image
                imfinfo( fileName );
               
                isImageFormat = true;
                
            catch exception %#ok<NASGU>
                % imfinfo does not recognize this file, don't error
                % since it is most likely a valid multimedia file
            end
            
            if isImageFormat
                % If imfinfo does not error, then show this error
                error(message('MATLAB:audiovideo:VideoReader:unsupportedImage'));
            end
        end
        
        function fileDesc = translateDescToLocale(fileExtension)
            switch upper(fileExtension)
                case 'M4V'
                    fileDesc = getString(message('MATLAB:audiovideo:VideoReader:formatM4V'));
                case 'MJ2'
                    fileDesc = getString(message('MATLAB:audiovideo:VideoReader:formatMJ2'));
                case 'MOV'
                    fileDesc = getString(message('MATLAB:audiovideo:VideoReader:formatMOV'));
                case 'MP4'
                    fileDesc = getString(message('MATLAB:audiovideo:VideoReader:formatMP4'));
                case 'MPG'
                    fileDesc = getString(message('MATLAB:audiovideo:VideoReader:formatMPG'));
                case 'OGV'
                    fileDesc = getString(message('MATLAB:audiovideo:VideoReader:formatOGV'));
                case 'WMV'
                    fileDesc = getString(message('MATLAB:audiovideo:VideoReader:formatWMV'));
                otherwise
                    % This includes formats such as AVI, ASF, ASX.
                    fileDesc = getString(message('MATLAB:audiovideo:VideoReader:formatGeneric', upper(fileExtension)));
            end
        end
    end
    
    %------------------------------------------------------------------
    % Helpers
    %------------------------------------------------------------------
    methods (Access='private', Hidden)

        function init(obj, fileName)
            
            % Properly initialize the object on construction or load.
            
            % Expand the path, using the matlab path if necessary
            fullName = audiovideo.internal.absolutePathForReading(...
                fileName, ...
                'MATLAB:audiovideo:VideoReader:fileNotFound', ...
                'MATLAB:audiovideo:VideoReader:FilePermissionDenied');

            VideoReader.errorIfImageFormat(fullName);
            
            % Save constructor arg for load.
            obj.ConstructorArgs = fullName;
            
            % Create underlying implementation.
            try
               obj.VideoReaderImpl = audiovideo.mmreader(fullName);
            catch exception
               VideoReader.handleImplException( exception );
            end
            
            if obj.EnableFrameCounting
                obj.populateNumFrames()
            
                % NumberOfFrames property is set to empty if it cannot be
                % determined by from the video. Generate a warning in this
                % case.
                if isempty(obj.NumberOfFrames)
                    warnState=warning('off','backtrace');
                    c = onCleanup(@()warning(warnState));
                    warning(message('MATLAB:audiovideo:VideoReader:unknownNumFrames'));
                end
            end
        end
        
        function impl = getImpl(obj)
            impl = obj.VideoReaderImpl;
        end
        
        function value = getImplValue(obj, propName)
            value = obj.getImpl().(propName);
        end
        
        
        function setImplValue(obj, propName, value) %#ok<INUSD>
            % All underlying properties are read only. Make the error 
            % the same as a standard MATLAB error when setting externally.
            % TODO: Remove when g449420 is done and used when calling 
            % set() in the constructor.
            err = MException('MATLAB:class:SetProhibited',...
                             'Setting the ''%s'' property of the ''%s'' class is not allowed.',...
                             propName, class(obj));
            throwAsCaller(err);
        end
        
        function [headings, indices] = getCategoryInfo(obj, propNames)
            % Returns headings and property indices for each category.
            headings = {'General Settings' 'Video Settings', 'Audio Settings'};
            indices = {[] [] []};
            for pi=1:length(propNames)
                propInfo = findprop(getImpl(obj), propNames{pi});
                if isempty(propInfo) || strcmpi(propInfo.Category, 'none')
                    category = 'general';
                else
                    category = propInfo.Category;
                end
                switch category
                    case 'general'
                        indices{1}(end+1) = pi;
                    case 'video'
                        indices{2}(end+1) = pi;
                    case 'audio'
                        indices{3}(end+1) = pi;
                end
            end
        end
    end
    
    methods (Hidden)
        function settableProps = getSettableProperties(obj)
            % Returns a list of publically settable properties.
            % TODO: Reduce to fields(set(obj)) when g449420 is done.
            settableProps = {};
            props = fieldnames(obj);
            for ii=1:length(props)
                p = findprop(obj, props{ii});
                if strcmpi(p.SetAccess,'public')
                    settableProps{end+1} = props{ii}; %#ok<AGROW>
                end
            end
        end
    end
        %}
end
