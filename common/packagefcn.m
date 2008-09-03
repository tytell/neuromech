function packagefcn(fcns,dest)

if (nargin == 0),
    [fn,pn] = uigetfile('*.m','Select function(s) to package', 'MultiSelect','on');
    
    if (~iscell(fn)),
        fn = {fn};
    end;
    
    fcns = cell(size(fn));
    for i = 1:length(fn),
        fcns{i} = fullfile(pn,fn{i});
    end;
    
    dest = uigetdir('','Select destination direction');
end;

package = {};
for i = 1:length(fcns),
    [pn,fn,ext] = fileparts(fcns{i});
    cd(pn);
    
    dep = depfun(fn);
    
    builtin = strmatch(matlabroot,dep);
    isbuiltin = false(size(dep));
    isbuiltin(builtin) = true;
    
    package = union(package,dep(~isbuiltin));
end;

for i = 1:length(package),
    copyfile(package{i},dest);
end;


