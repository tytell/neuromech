function sel = selectgroups(varargin)
% function sel = selectgroups(gp1,gp2,gp3...)
%   or           selectgroups(groups)
%
% Selects unique combinations of different groups using a GUI.  gp1 etc are
% variables that have the same size.  selectgroups attempts to extract the
% name of the group from the variable name, but it doesn't always work.
% groups is either a cell array of variables or a structure.  All elements
% of the cell array or structure need to have the same size.
%
% Can also pass an option, 'groupnames', followed by a cell string, to
% specify the names of the groups.

gpnames = {};
gap = 0.02;
if (isstruct(varargin{1})),
    S = varargin{1};
    gpnames = fieldnames(S);
    groups = struct2cell(S);
    p = 2;
elseif (iscell(varargin{1})),
    groups = varargin{1};
    p = 2;
else
    p = 1;
    
    nd = ndims(varargin{1});
    sz = size(varargin{1});
    
    while ((p <= nargin) && (ndims(varargin{p}) == nd) && ...
            all(size(varargin{p}) == sz)),
        gpnames{p} = inputname(p);
        p = p+1;
    end;
end;

while (p <= nargin)
    switch lower(varargin{p}),
        case 'groupnames',
            gpnames = varargin{p+1};
            p = p+2;
            
        otherwise,
            error('Unrecognized option %s', varargin{p});
    end;
end;

%get rid of empty groups
good = ~cellfun(@isempty,groups);
groups = groups(good);
if (~isempty(gpnames)),
    gpnames = gpnames(good);
end;

ngp = length(groups);
nd = ndims(groups{1});
sz = size(groups{1});

%open the figure and set up the data structure
fig = openfig(mfilename, 'new');
data = guihandles(fig);
pos0 = get(data.ListBox, 'Position');
w = (pos0(3)-(ngp-1)*gap)/ngp;
gpind = zeros(numel(groups{1}),ngp);

nlevels = zeros(1,ngp);
hlist = zeros(1,ngp);
selvals = cell(1,ngp);
for i = 1:ngp,
    if ((ndims(groups{i}) ~= nd) || any(size(groups{i}) ~= sz)),
        error('All groups must be the same size');
    end;

    %get the unique, non-nan values from each group
    good = ~isnan(groups{i}(:));
    if (any(good)),
        [gpvals,q,gpind(good,i)] = unique(groups{i}(good));
        
        %leave one NaN value at the end, if there are any
        if (any(~good)),
            gpvals(end+1) = NaN;
            gpind(~good,i) = length(gpvals);
        end;
    else
        gpvals = NaN;
        gpind(:,i) = 1;
    end;
    nlevels(i) = length(gpvals);
    
    %add a listbox, based on duplicating one initial prototype
    hlist(i) = copyobj(data.ListBox,fig);
    
    %set up the strings to go in the listbox
    str = cell(nlevels(i),1);
    if (isnumeric(gpvals)),
        for j = 1:nlevels(i),
            str{j} = num2str(gpvals(j));
        end;
    elseif (iscellstr(gpvals)),
        str{i} = gpvals;
    elseif (ischar(gpvals)),
        for j = 1:nlevels(i),
            str{j} = gpvals(j);
        end;
    elseif (islogical(gpvals)),
        if (nlevels(i) == 1) 
            if (gpvals),
                str{1} = 'true';
            else
                str{1} = 'false';
            end;
        else
            str{1} = 'false';
            str{2} = 'true';
        end;
    end;
    
    %and set the position
    l = (i-1)*(w + gap);
    pos1 = get(hlist(i),'Position');
    pos1(1) = pos1(1) + l;
    pos1(3) = w;
    
    %and the initial selection
    if (nlevels(i) == 1),
        mx = 1;
        val = 1;
        selvals{i} = 1;
    else
        mx = nlevels(i)+1;
        val = [];
        selvals{i} = [];
    end;
    set(hlist(i), 'Position',pos1, 'String',str, ...
        'Max',mx,'Min',1, 'Value',val, 'Callback',{@sgListClick,i});
    
    if (~isempty(gpnames) && ~isempty(gpnames{i})),
        hlabel = copyobj(data.GroupLabel,fig);
    
        pos1 = get(hlabel,'Position');
        pos1(1) = l;
        pos1(3) = w;
        set(hlabel, 'Position',pos1, 'String',gpnames{i});
    end;
    
    hnum(i) = copyobj(data.NumLabel,fig);
    
    pos1 = get(hnum(i),'Position');
    pos1(1) = l;
    pos1(3) = w;
    if (nlevels(i) == 1),
        str = num2str(sum(gpind(:,i) == 1));
    else
        str = '0';
    end;
    set(hnum(i), 'Position',pos1, 'String',str);
end;

%get rid of the prototype listbox and label
delete(data.ListBox);
delete(data.GroupLabel);

%and save the ones we just created
data.hlist = hlist;
data.hnum = hnum;
data.gpind = gpind;
data.selvals = selvals;

%set up a status display
set(data.TotalNum,'String',sprintf('Total number: %d',size(gpind,1)));
set(data.SelectedNum,'String',sprintf('Selected: %d',0));

%and appropriate callbacks
set(data.OKButton,'Callback',@sgOKCancel);
set(data.CancelButton,'Callback',@sgOKCancel);
guidata(fig,data);

%and run
uiwait(fig);
if (ishandle(fig)),
    sel = get(fig,'UserData');
    if (~isempty(sel)),
        sel = reshape(sel,sz);
    end;
    delete(fig);
else
    sel = [];
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nsel,nselindiv,sel] = sgGetNumSel(gpind,selvals)
%set up a logical matrix that defines which items are selected

%for each group, and with the selected values
sel = true(size(gpind,1),1);
nselindiv = zeros(1,length(selvals));
for i = 1:size(gpind,2),
    sel1 = false(size(gpind,1),1);
    for j = 1:length(selvals{i}),
        sel1 = sel1 | (gpind(:,i) == selvals{i}(j));
    end;
    nselindiv(i) = sum(sel1);
    sel = sel & sel1;
end;
nsel = sum(sel);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sgListClick(obj, eventdata, gpnum)
%handle a click on the list

data = guidata(obj);
val = get(obj, 'Value');
data.selvals{gpnum} = val;
[nsel,nselindiv,sel] = sgGetNumSel(data.gpind,data.selvals);
set(data.SelectedNum,'String',sprintf('Selected: %d',nsel));
set(data.Figure,'UserData', sel);
for i = 1:length(nselindiv),
    set(data.hnum(i),'String',num2str(nselindiv(i)));
end;
guidata(obj,data);

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sgOKCancel(obj, eventdata)
data = guidata(obj);
if (obj == data.OKButton),
    [nsel,nselindiv,sel] = sgGetNumSel(data.gpind,data.selvals);
    set(data.SelectedNum,'String',sprintf('Selected: %d',nsel));
    set(data.Figure,'UserData', sel);
else
    set(data.Figure, 'UserData', []);
end;
uiresume;            
    
   