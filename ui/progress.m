function progress(i,n, msg, varargin)

opt.show = 'percent';
opt.percent = 0.1;
opt.time = 10;

opt = parsevarargin(opt,varargin, 4);

global prg_ticstart;
global prg_numsteps;
global prev_elap;

if ((i == 0) || (i == 1))
    prg_ticstart = tic;
    fprintf('%s\n', msg);
    prg_numsteps = n;
    prev_elap = 0;
elseif (i <= prg_numsteps)
    elap = toc(prg_ticstart);
    switch opt.show
        case 'all'
            show = true;
        case 'percent'
            pctsteps = opt.percent*prg_numsteps;
            show = floor((i-1)/pctsteps) < floor(i/pctsteps);
        case 'time'
            show = elap - prev_elap > opt.time;
    end
    
    if (show)
        remain = (elap/i) * (prg_numsteps-i);
        fprintf('%d/%d steps (%d%%). %s elapsed; %s remaining\n', ...
            i,prg_numsteps, round((i/prg_numsteps)*100), ...
            formatseconds(elap), formatseconds(remain));
    end
end
            
    