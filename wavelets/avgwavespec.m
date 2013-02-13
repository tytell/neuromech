function Wavg = avgwavespec(scales,W,smoothscale,smoothsamples)

if (all(diff(scales) < 0)),
    W = flipud(W);
    scales = scales(end:-1:1);
    flip = true;
else
    flip = false;
end;

a = scales;
Nscales = length(a);

if ((smoothscale > 0) && (length(a) > 1))
    octaves = log2(a);
    
    Nvoiceall = 1./abs(diff(octaves));
    Nvoice = mean(Nvoiceall);
    
    wa = round(smoothscale * Nvoice);
    wa2 = 2*wa+1;

    Wavg = NaN(size(W));
    
    scaleavg = filter(ones(wa2,1)/wa2,1, W);
    Wavg(wa+1:end-wa, :,:) = scaleavg(wa2:end, :,:);
else
    Wavg = W;
    wa = 0;
end;

wb = round(smoothsamples * a);
wb2 = 2*wb+1;
for i = wa+1:Nscales-wa,
    if (wb2(i) > 1),
        timeavg = filter(ones(1,wb2(i))/wb2(i),1, Wavg(i,:,:), [], 2);
        Wavg(i,wb(i)+1:end-wb(i),:) = timeavg(:,wb2(i):end,:);
    end;
end;
    
if (flip)
    Wavg = flipud(Wavg);
end;

        
    
    

