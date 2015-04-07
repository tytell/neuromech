function progress(i,n, msg, varargin)

opt.show = 'percent';
opt.percent = 0.1;
opt.time = 10;

opt = parsevarargin(opt,varargin, 4);

global prg_ticstart;
global prg_totalsteps;
global prg_numsteps;
global prev_elap;

if (length(i) > 1) 
    if ((nargin >= 3) && (length(n) == length(i)))
        prg_numsteps = n;
        sz = [1; cumprod(n(1:end-1))];
    else
        sz = [1; cumprod(prg_numsteps(1:end-1))];
    end
    if all(i == 0)
        ind = 0;
    elseif all((i >= 1) & (i <= prg_numsteps))
        ind = sum((i(:)-1).*sz) + 1;
    end
    i = ind;
end

if (i == 0) || ((i == 1) && (nargin >= 3))
    prg_ticstart = tic;
    fprintf('%s\n', msg);
    prg_totalsteps = prod(n);
    prev_elap = 0;
elseif (i <= prg_totalsteps)
    elap = toc(prg_ticstart);
    switch opt.show
        case 'all'
            show = true;
        case 'percent'
            pctsteps = opt.percent*prg_totalsteps;
            show = floor((i-1)/pctsteps) < floor(i/pctsteps);
        case 'time'
            show = elap - prev_elap > opt.time;
    end
    
    if (show || (i <= 1))
        remain = (elap/i) * (prg_totalsteps-i);
        fprintf('%d/%d steps (%d%%). %s elapsed; %s remaining\n', ...
            i,prg_totalsteps, round((i/prg_totalsteps)*100), ...
            formatseconds(elap), formatseconds(remain));
    end
end
            
    