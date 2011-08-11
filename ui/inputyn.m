function tf = inputyn(str,varargin)
%INPUTYN  Asks a yes or no question and returns true or false.
%
% function tf = inputyn(str,default)
% Returns true if the user types anything with a 'y' or 'Y' in it
% Returns default if the user doesn't type anything.
% Default is optional and defaults to true.
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell

opt.default = true;
opt.strict = false;
opt.showdefault = true;
opt.matchallcallers = false;
opt.quiet = false;
opt.clearsaved = false;
opt.saveanswer = {};
opt.requireanswer = false;

helpstr = {'''y'' or ''yes'' - yes', ...
    '''n'' or ''no'' - no', ...
    '''all yes'' or ''yes all'' - always yes to this question', ...
    '''all no'' or ''no all'' - always no to this question', ...
    '''clear'' - clear any saved answers', ...
    '''h'' or ''help'' or ''?'' - display this text'};

persistent savedanswers;

opt = parsevarargin(opt,varargin,2, 'typecheck',false);

if (isempty(savedanswers))
    savedanswers = {};
end;

if (opt.clearsaved)
    if (ischar(opt.clearsaved))
        keep = ~strcmp(str,savedanswers(:,1));
        savedanswers = savedanswers(keep,:);
    else
        savedanswers = {};
    end;
    return;
elseif (~isempty(opt.saveanswer))
    %programmatically save the answer to a question
    savedanswers(end+(1:size(opt.saveanswer,1)),1:3) = opt.saveanswer;
    return;
end;
    
if (opt.showdefault)
    if (opt.strict)
        showstr = [str ' (type ''yes'' to confirm) '];
    elseif (opt.default)
        showstr = [str ' (y,n,all yes,all no,clear,?) (default = ''y'') '];
    else
        showstr = [str ' (y,n,all yes,all no,clear,?) (default = ''n'') '];
    end;
end;

issaved = false;
st = dbstack;
if (length(st) >= 2)
    callfunc = st(2).name;
else
    callfunc = '';
end;

if (~opt.strict && ~isempty(savedanswers))
    ind = find(strcmp(str,savedanswers(:,1)));
    if (length(ind) == 1)
        if (opt.matchallcallers || strcmp(savedanswers{ind,2},callfunc))
            tf = savedanswers{ind,3};
            if (tf)
                showans = 'yes';
            else
                showans = 'no';
            end;
            
            fprintf('%s Using saved answer = %s\n', str, showans);
            issaved = true;
        end;
    end;
end;

if (~issaved)
    done = false;
    
    while (~done)
        yn = input(showstr,'s');
        if (opt.strict)
            tf = strcmpi(yn,'yes');
            done = true;
        else
            switch lower(yn)
                case {'y','yes'}
                    tf = true;
                    done = true;
                    
                case {'n','no'}
                    tf = false;
                    done = true;
                    
                case {'all yes', 'yes all'},
                    savedanswers(end+1,1:3) = {str,callfunc,true};
                    tf = true;
                    done = true;
                    
                case {'all no','no all'}
                    savedanswers(end+1,1:3) = {str,callfunc,false};
                    tf = false;
                    done = true;
                    
                case 'clear'
                    fprintf('Cleared old saved answers.\n');
                    savedanswers = {};
                    done = false;
                    
                case {'h','help','?'}
                    fprintf('%s\n', helpstr{:});
                    done = false;
                    
                otherwise
                    if (~opt.requireanswer)
                        tf = opt.default;
                        done = true;
                    end;
            end;
        end;
    end;
end;
