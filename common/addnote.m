function addnote(aStr, noteNum)
% addnote(aStr)
% adds a note the cell array workspaceNotes.
%
% addnote
% displays the notes in the array workspaceNotes.
%
% addnote([])
% clears the variable 'workspaceNotes' in the base workspace (or function
% workspace if addnote is called from within a function).
%
% addnote([], noteNum)
% deletes the note noteNum
%
% addnote(aStr, noteNum)
% overwrite the note noteNum
%
% This function will add a note to a function workspace as well
% if it is called from within a function. 
%
% Peter Madden, July 6, 2005.
%

ws = evalin('caller', 'myWorkspace', '[]');

if isfield(ws, 'Notes')
    wsNotes = ws.Notes;
else
    wsNotes = [];
end

if nargin == 0    
    if isempty(wsNotes)
        disp('No notes saved in myWorkspace.');
    else
        for i = 1:length(wsNotes)
            disp(wsNotes{i});
        end
    end
elseif nargin == 1      % add a string or erase all if aStr = []
    if ~isempty(aStr)   % add a string
        dStr = datestr(now, 31);
        wsStr = [dStr, '  : ', aStr, ' (',num2str(length(wsNotes)+1),')'];
        if isempty(wsNotes)
            wsNotes{1} = wsStr;
        else
            wsNotes{end + 1} = wsStr;
        end
        evalin('caller', 'if ~exist(''myWorkspace''), myWorkspace.Notes = {}; end');
        assignin('caller', 'myWorkspaceTempVar', wsNotes);
        evalin('caller', 'myWorkspace.Notes = myWorkspaceTempVar;');
        evalin('caller', 'clear myWorkspaceTempVar;');
    elseif isempty(aStr)    % clear all Notes
        button = questdlg('Clear the workspaceNotes variable?',...
                    'Clear workspaceNotes?','Yes','No','Cancel','No');
        if strcmp(button,'Yes')
            evalin('caller', 'if ~exist(''myWorkspace''), myWorkspace = []; end');
            evalin('caller', 'myWorkspace.Notes=[];');
        end    
    end
elseif nargin == 2  % overwrite a string or delete a specific string
    if ~isempty(aStr)   % --- overwrite --- 
        dStr = datestr(now, 31);
        wsStr = [dStr, '  : ', aStr, ' (',num2str(noteNum),')'];        
        wsNotes{noteNum} = wsStr;
        assignin('caller', 'myWorkspaceTempVar', wsNotes);                
        evalin('caller', 'myWorkspace.Notes = myWorkspaceTempVar;');
        evalin('caller', 'clear myWorkspaceTempVar;');        
        
    elseif isempty(aStr)    % --- delete ---
        for i = noteNum +1:length(wsNotes)
            wsNotes(i-1:end-1) = wsNotes(i:end);
        end
        wsNotes = wsNotes(1:end-1);
        evalin('caller', 'if ~exist(''myWorkspace''), myWorkspace.Notes = {}; end');   
        assignin('caller', 'myWorkspaceTempVar', wsNotes);        
        evalin('caller', 'myWorkspace.Notes = myWorkspaceTempVar;');
        evalin('caller', 'clear myWorkspaceTempVar;');
    end
end
        
