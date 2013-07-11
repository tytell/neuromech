function colorby(h,gp, varargin)

opt.color = 'bgrmcyk';
opt = parsevarargin(opt,varargin, 3);

[gpnm,~,gp] = unique(gp,'stable');

for i = 1:length(gpnm)
    isgp = gp == i;
    
    j = mod(i-1,length(opt.color)) + 1;
    set(h(isgp),'Color',opt.color(j));
end
