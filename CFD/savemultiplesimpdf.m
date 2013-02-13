function savemultiplesimpdf(outfile,datafiles, varargin)
% function savemultiplesimpdf(outfile,datafiles, varargin)
% Produces a large PDF file containing summary data from many simulations.

opt.showfig = 1:8;

opt = parsevarargin(opt,varargin,2);

nfig = length(opt.showfig);
fig = -1*ones(nfig,length(datafiles));

for f = 1:length(datafiles),
    plotsimdata(datafiles{f},'showfig',opt.showfig);

    for i = 1:nfig,
        fig(i,f) = copyobj(figure(opt.showfig(i)),0);
    end;
end;

delete(outfile);
for i = 1:nfig,
    for f = 1:length(datafiles),
        export_fig(outfile, '-pdf','-append',fig(i,f));
        close(fig(i,f));
    end;
end;

    