function [P,F,df1,df2] = angequalmeans(ang1,ang2,method)
% ANGEQUALMEANS - Test for equal means in two groups of angles
%
%   [P,F,df1,df2] = angequalmeans(ang1,ang2,method)
%
% Test for equal means in the two groups ang1 and ang2 using one of
% the following methods:
%     'watson-williams','ww','parametric': Parametric test for equality of means
%     'watson','w','nonparametric': Watson's U^2 test, a nonparametric test
%              for a common distribution
%
% NB: anganova is more powerful and will perform two sample tests (although
% not the Watson U^2)
%
% SEE ALSO
%   ANGANOVA, ANGMEAN

ang1 = shiftdim(ang1(isfinite(ang1)));
ang2 = shiftdim(ang2(isfinite(ang2)));

switch lower(method),
 case {'watson-williams','ww','parametric'},
     %from Fisher, 1993 p.126
     %(assumes equal concentrations, all >= 2)
     angx1 = cos(ang1);
     angy1 = sin(ang1);
     angx2 = cos(ang2);
     angy2 = sin(ang2);

     n1 = length(ang1);
     n2 = length(ang2);
     N = n1 + n2;

     angx = mean([angx1; angx2]);
     angy = mean([angy1; angy2]);
     angx1 = mean(angx1);
     angy1 = mean(angy1);
     angx2 = mean(angx2);
     angy2 = mean(angy2);

     R1 = n1*sqrt(angx1^2 + angy1^2);
     R2 = n2*sqrt(angx2^2 + angy2^2);
     R = N*sqrt(angx^2 + angy^2);

     F = (N-2)*(R1+R2-R)/(N-(R1+R2));

     df1 = 1;
     df2 = N-2;
     P = 1 - fcdf(F,df1,df2);

    case {'watson','w','nonparametric'},
        %from Batschelet, 1981, p.114
        %tests differences in mean or concentration simultaneously
        n1 = length(ang1);
        n2 = length(ang2);
        N = n1 + n2;

        [ang,ind] = sort(cat(1,ang1,ang2));

        rank1 = [ones(1,n1) zeros(1,n2)];
        rank2 = [zeros(1,n1) ones(1,n2)];

        rank1 = cumsum(rank1(ind))/n1;
        rank2 = cumsum(rank2(ind))/n2;

        d = rank1 - rank2;

        U2 = n1*n2/N^2*(sum(d.^2) - sum(d)^2/N);

        P = watsonU2(U2,n1,n2);
        F = U2;
        df1 = n1;
        df2 = n2;
end;

  
  


