function [opt,rest,indrest] = parsevarargin(opt,varargs,varargin)
%PARSEVARARGIN    Simple yet powerful way of handling optional arguments
% opt = parsevarargin(opt,varargs,[firstoptnum],...)
%  or
% [opt,rest] = parsevarargin(...)
%
% Parses varargin for another function and fills the opt structure.
% Provides a simpler structure than the inputParser class, but still has
% many of the same options.  opt is a structure with field names that are
% the option names and values that are the defaults for each option.
% varargs is the varargin from the calling function.  firstoptnum is an
% optional scalar, defining the index of the first optional argument in the
% calling function, so that parsevarargin can return errors that correctly
% reflect the number of required and optional arguments to the calling
% function.  This equivalent to passing the 'firstoptionnumber'.
%
% Usually tries to match short names, provided they are unambiguous (ie, 'err'
% matches 'error', provided there's not another option called
% 'errorfunc').
%
% Example:
%    function testfunc(a,b,varargin)
%
%    opt.debug = false;
%    opt.param = 12;
%    opt.method = 'eat';
%
%    opt = parsevarargin(opt,varargin, 3, ...
%             'multival',{'method',{'eat','chew'}}, ...
%             'synonyms',{{'debug',{'dbg'}}});
%
%  So all of the following would be valid calls to testfunc:
%    testfunc(a,b)
%       -> opt is set to default values
%    testfunc(a,b,'debug')
%       -> opt.debug = true
%    testfunc(a,b,'dbg',true)
%       -> opt.debug = true
%    testfunc(a,b,'method','chew')
%       -> opt.method = 'chew'
%    testfunc(a,b,'chew')
%       -> opt.method = 'chew'
%    testfunc(a,b,'p',20)
%       -> opt.param = 20
%
% Options
%
%   'allowno' - Allows "no" options.  For example, if there was an option 'showerror' that
%   was normally true, then if 'allowno' was true, the option 'noshowerror' would also be
%   valid.
%
%   'typecheck' - Checks to make sure the arguments for each option match the type of the
%   default value.
%
%   'firstoptionnumber' - Gives the parameter number of the first option, so that
%   parameter numbers can be returned correctly.  For example, to parse the options for
%   curvature(x,y,varargin), use parsevarargin(...,'firstoptionnumber',3).
%
%   'multival' - Handles multiple valued options.  The argument is a cell array in which
%   the first value is the name of the option and the second value is a cell string of
%   potential values for that option.  Multiple rows can be used for different multiple
%   valued options.  For example, in curvature.m, the option 'method' can be 'spline' or
%   'discrete'.  With parsevarargin(...,'multival',{'method',{'spline','discrete'}}) the
%   user could write curvature(...,'spline').
%
%   'allowunknown' - Allows unknown options.  For options that aren't listed in opt, it
%   won't throw an error, it will just add the option and value to the opt structure.
%   Unknown options must have an option name and an argument.
%
%   'leaveunknown' - Parses the options that match and returns the rest as
%   a second output.
%
%   'synonyms' - Allows a list of synonynms for argument names, as follows:
%   {{'option',{'synonym1','synonym2'}}, {'argument',{'arg','synonym3'}}.
%   In this case, the user could write 'argument',3 or 'arg',3 and they
%   would both map to opt.argument.
%
%   'exact' - Force only exact matches to the option names.  
%
% See also INPUTPARSER.

% Mercurial revision hash: $Revision$ $Date$
% See https://bitbucket.org/tytell/matlab/wiki/Home
% Copyright (c) 2011, Eric Tytell <tytell at jhu dot edu>


%our own defaults
pvopt.allowno = false;
pvopt.typecheck = true;
pvopt.firstoptionnumber = [];
pvopt.multival = [];
pvopt.allowunknown = false;
pvopt.leaveunknown = false;
pvopt.synonyms = {};
pvopt.exact = true;
firstoptnum = [];

%process those options
if ((nargin >= 3) && isnumeric(varargin{1}) && (numel(varargin{1}) == 1)),
    firstoptnum = varargin{1};
    p = 2;
else
    p = 1;
end;

while (p <= length(varargin)),
    switch lower(varargin{p}),
      case {'firstoptionnumber','multival','synonyms'},
        pvopt.(lower(varargin{p})) = varargin{p+1};
        p = p+2;
        
      case {'typecheck','allowno','allowunknown','leaveunknown','exact'},
        if ((p+1 <= length(varargin)) && islogical(varargin{p+1})),
            pvopt.(lower(varargin{p})) = varargin{p+1};
            p = p+2;
        else
            pvopt.(lower(varargin{p})) = true;
            p = p+1;
        end;
        
      otherwise,
        error('pvopt:UnrecognizedOptionToParsevarargin','Unrecognized option %s',varargin{p});
    end;
end;

if (~isempty(pvopt.firstoptionnumber) && isempty(firstoptnum))
    firstoptnum = pvopt.firstoptionnumber;
end;
if (isempty(firstoptnum)),
    firstoptnum = 1;
end;

%get the option names
optnames = fieldnames(opt);
%and their default values
defaults = struct2cell(opt);

%check for logical options (those that have true or false as their defaults)
logicalopt = cellfun(@islogical,defaults);
%construct the "no" options
if (pvopt.allowno),
    noopt = cellfun(@(x) (strcat('no',x)), optnames(logicalopt),'UniformOutput',false);
else
    noopt = {};
end;

%set up synonyms
synnames = [];
synoptind = [];       % the option name that the synonym refers to
for i = 1:length(pvopt.synonyms)
    if (~ischar(pvopt.synonyms{i}{1}) || ~iscellstr(pvopt.synonyms{i}{2}))
        error('pvopt:BadSynonymSpec','Synonym array is not structured right');
    end;
    ind1 = find(strcmpi(pvopt.synonyms{i}{1},optnames));
    if (length(ind1) ~= 1)
        error('pvopt:BadSynonymSpec','Couldn''t find base name for synonym %s', ...
            pvopt.synonyms{i}{1});
    end;
    synnames1 = pvopt.synonyms{i}{2};
    
    synnames = [synnames synnames1];        %#ok
    
    synoptind = [synoptind ind1*ones(1,length(synnames1))];    %#ok
end;

isrest = false(size(varargs));

%run through the varargs passed in
p = 1;
while (p <= length(varargs)),
    ismatch = false;
    
    %check for a char argument
    if (ischar(varargs{p}) && (size(varargs{p},1) == 1))
        if (pvopt.exact)
            match = find(strcmpi(varargs{p},optnames));
            syn = find(strcmpi(varargs{p},synnames));
        else
            match = find(strncmpi(varargs{p},optnames,length(varargs{p})));
            syn = find(strncmpi(varargs{p},synnames,length(varargs{p})));
        end;
        
        if (length(syn) == 1)
            match = synoptind(syn);
        end;
        if ((length(match) == 1) && ~logicalopt(match))
            %handle normal options
            
            %check to make sure there's an argument
            if (p+1 > length(varargs)),
                ME = MException('pvopt:missingoptionarg','Option %s is missing an argument',...
                    lower(varargs{p}));
                throwAsCaller(ME);
            elseif (pvopt.typecheck && ~strcmp(class(varargs{p+1}), ...
                    class(opt.(optnames{match})))),
                %do the type check, if necessary
                ME = MException('pvopt:optiontypemismatch',...
                    'Option %s requires an argument of type %s',...
                    lower(varargs{p}),class(opt.(optnames{match})));
                throwAsCaller(ME);
            else
                %save the option
                opt.(optnames{match}) = varargs{p+1};
                ismatch = true;
            end;
            p = p+2;
        elseif ((length(match) == 1) && logicalopt(match))
            %handle logical options
            
            %check for a true/false argument
            if ((p+1 <= length(varargs)) && islogical(varargs{p+1})),
                opt.(optnames{match}) = varargs{p+1};
                p = p+2;
            else
                opt.(optnames{match}) = true;
                p = p+1;
            end;
            ismatch = true;
        elseif (isempty(match) && ~isempty(pvopt.multival)),
            %handle multiple valued options
            for i = 1:size(pvopt.multival,1),
                if (ismember(lower(varargs{p}),pvopt.multival{i,2})),
                    opt.(pvopt.multival{i,1}) = lower(varargs{p});
                    ismatch = true;
                    p = p+1;
                    break;
                end;
            end;
        elseif (pvopt.allowno),
            %handle "no" options
            noval = varargs{p};
            if (pvopt.exact)
                matchno = strcmpi(noval,noopt);
            else
                matchno = strncmpi(noval,noopt,length(noval));
            end;
            if (sum(matchno) == 1)
                optname = noopt{matchno};
                optname = optname(3:end);
                opt.(optname) = false;
                ismatch = true;
                p = p+1;
            end;
        end;
    end;
    
    %if we don't have a match, save the unmatched option or throw an error
    if (~ismatch)
        if (pvopt.leaveunknown),
            isrest(p) = true;
            p = p+1;
        elseif (pvopt.allowunknown),
            %unmatched options must have arguments
            if (p+1 > length(varargs)),
                ME = MException('pvopt:missingoptionarg','Option %s is missing an argument',...
                                varargs{p});
                throwAsCaller(ME);
            else
                opt.(varargs{p}) = varargs{p+1};
                p = p+2;
            end;
        else
            if (ischar(varargs{p})),
                ME = MException('pvopt:unrecognizedoption',...
                    'Unrecognized option %s at parameter %d',...
                    lower(varargs{p}), p+firstoptnum-1);
            else
                ME = MException('pvopt:unrecognizedoption',...
                    'Unrecognized option at parameter %d',...
                    p+firstoptnum-1);
            end;
            throwAsCaller(ME);
        end;
    end;
end;

if (pvopt.leaveunknown),
    rest = varargs(isrest);
    indrest = find(isrest);
else
    rest = {};
end;

        
            
            
        
        
            
