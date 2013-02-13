function varargout = findactivationwave(actl,actr)
% function indact = findactivationwave(actl,actr)
% Returns an indact structure, like indpeak from analyzekinematics, which
% tracks the position of the activation waves along the body.  indact is
% a matrix with the same number of rows as points on the body, and the same
% number of columns as activation waves.  The values in the matrix are the
% frame number when the beginning of the activation wave was at the
% corresponding point (=row) in the matrix.

nfr = size(actl,2);
npt = size(actl,1);

%outline connected regions in the time-space activation matrix
lwave = bwboundaries(actl,'noholes');
indactl = NaN(npt,length(lwave));
indactoffl = NaN(npt,length(lwave));
for i = 1:length(lwave),
    %find the minimum frame at which a specific point is active
    actlon1 = accumarray(lwave{i}(:,1),lwave{i}(:,2), [npt 1], @min, NaN);
    %and the maximum frame
    actloff1 = accumarray(lwave{i}(:,1),lwave{i}(:,2), [npt 1], @max, NaN);
    
    %put those in out indact matrices
    indactl(:,i) = actlon1;
    indactoffl(:,i) = actloff1;
end;

%same thing for right side activation
rwave = bwboundaries(actr,'noholes');
indactr = NaN(npt,length(rwave));
indactoffr = NaN(npt,length(rwave));
for i = 1:length(rwave),
    actron1 = accumarray(rwave{i}(:,1),rwave{i}(:,2), [npt 1], @min, NaN);
    actroff1 = accumarray(rwave{i}(:,1),rwave{i}(:,2), [npt 1], @max, NaN);

    indactr(:,i) = actron1;
    indactoffr(:,i) = actroff1;
end;

if (nargout == 4),
    varargout = {indactl,indactoffl, indactr,indactoffr};
else
   %put the two matrices together, assuming they alternate
   indact = NaN(npt,size(indactl,2)+size(indactr,2));
   indactoff = NaN(npt,size(indactl,2)+size(indactr,2));

   %which side comes first?
   if (max(indactl(:,1)) < max(indactr(:,1))),
       l1 = 1;
       r1 = 2;
   else
       l1 = 2;
       r1 = 1;
   end;
   indact(:,l1:2:end) = indactl;
   indact(:,r1:2:end) = indactr;
   indactoff(:,l1:2:end) = indactoffl;
   indactoff(:,r1:2:end) = indactoffr;

   isleft = false(1,size(indact,2));
   isleft(l1:2:end) = true;
   
   switch nargout,
     case 3,
       varargout = {indact,indactoff, isleft};
     case 2,
       varargout = {indact,indactoff};
     case 1,
       varargout = {indact};
       
     otherwise,
       error('Wrong number of outputs');
   end;
end;
