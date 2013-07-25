function colorby(h,gp, varargin)

opt.color = 'bgrmcyk';
opt = parsevarargin(opt,varargin, 3, 'typecheck',false);

[gpnm,~,gp] = unique(gp,'stable');

col = opt.color;
if (ischar(col))
    col = col(:);
end;

for i = 1:length(gpnm)
    isgp = gp == i;
    
    j = mod(i-1,size(col,1)) + 1;
    set(h(isgp),'Color',col(j,:));
end
