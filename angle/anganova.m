function varargout = anganova(ang, varargin)
% function P = anganova(ang, method, group)
% One way "ANOVA" type test on circular data.
% Available methods:
%     'watson-williams','ww': Parametric test for equality of means
%     'likelihood ratio','likelihood','lr': Parametric likelihood ratio
%          test for equality of means
%     'concentration','conc': Parametric test for equality of concentration
%     'mardia-watson-wheeler','mww','npdist': Nonparametric test for
%          equality of distribution (i.e., equal mean and concentration)
%     'median': Nonparametric test for equality of medians
%     'watson','npmean': Nonparametric test for equality of means
% If no output is requested, prints tables of test statistics and performs
% a Bonferroni-adjusted post hoc pairwise test afterwards, using the same
% method for the multisample test.
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell

alpha = 0.05;
group = [];
method = '';
quiet = false;
display = 'short';

i = 1;
while (i <= length(varargin)),
    if (isnumeric(varargin{i})),
        if ((ndims(varargin{i}) == ndims(ang)) && ...
                all(size(varargin{i}) == size(ang))),
            group = varargin{i};
            i = i+1;
        else
            error('anganova:badGroupParameter','Group parameter must be the same size as ang');
        end;
    elseif (ischar(varargin{i})),
        switch lower(varargin{i}),
            case 'alpha',
                alpha = varargin{i+1};
                i = i+2;
                
            case {'watson-williams','ww','likelihood ratio','likelihood','lr',...
                    'concentration','conc','mardia-watson-wheeler','mww','npdist',...
                    'median','watson','npmean'},
                method = varargin{i};
                i = i+1;
                
            case 'method',
                method = varargin{i+1};
                i = i+2;
                
            case 'quiet',
                quiet = true;
                i = i+1;
                
            case 'display',
                if ((i < length(varargin)) && ischar(varargin{i+1}) && ...
                    ismember(lower(varargin{i+1}),{'short','long','off'})),
                    display = lower(varargin{i+1});
                    i = i+2;
                else
                    quiet = false;
                    display = 'short';
                end;
                
            otherwise,
                error('anganova:badOption','Unrecognized option %s.',varargin{i});
        end;
    end;
end;

if (strcmpi(display,'off')),
    quiet = true;
end;

if (~isempty(group)),
    %rearrange ang so that each column is a group
    
    %eliminate NaN groups
    if (any(~isfinite(group(:)))),
        nonan = find(isfinite(group));
        ang = ang(nonan);
        group = group(nonan);
    end;

    %first set it up so that all ang0 > 0
    mm = min(ang(:));
    ang0 = ang - mm + 1;
    ang0(isnan(ang)) = 0;

    %now find the unique groups
    [gpind,gp] = grp2idx(group);
    gp = gp';

    N = accumarray(gpind,1);
    
    ang = NaN(N,length(gp));
    
    %run through all groups, putting each one in a column
    for i = 1:length(gp),
        isgrp = gpind == i;
        ang(1:sum(isgrp),i) = ang0(isgrp);
    end;

    %zeros will be at the end of the columns, so we set zeros to NaN,
    %remembering that we cleverly added on mm to ang0 so that all angles
    %were > 0
    ang(ang == 0) = NaN;
    %and now we subtract off mm again to return the angles back to where
    %they were
    ang = ang + mm - 1;
else
    gp = 1:size(ang,2);
end;

%check to make sure we have phase mod 2pi, not mod 1
if (range(ang) <= 1),
    warning('anganova:PhaseMod1','Phase appears to be mod 1.  Is this correct?');
end;

%get the average angle in each group
C = nansum(cos(ang));
S = nansum(sin(ang));
angmean = atan2(S,C);

%and the overall average angle
Call = sum(C);
Sall = sum(S);
angmeanall = atan2(Sall,Call);

%the R values
R = sqrt(C.^2 + S.^2);
Rall = sqrt(Call^2 + Sall^2);

%number of elements
ngp = size(ang,2);
n = sum(isfinite(ang));
nall = sum(isfinite(ang(:)));

%and the dispersion
k = angkappa(ang);
kall = angkappa(ang(:));

%print it all out in a nice table
if (~quiet && strcmp(display,'long')),
    fprintf('%8s %8s %8s %8s %8s %8s %8s %8s\n', 'Sample','C','S','R','N',...
        'mean','kappa','median');
    if (iscellstr(gp)),
        gptplt = '%8s';
    elseif (all(mod(gp,1) == 0)),
        gptplt = '%8d';
    else
        gptplt = '%8.3f';
    end;
    fprintf([gptplt ' %8.3f %8.3f %8.3f %8d %8.1f %8.2f %8.1f\n'],...
        [gp; C; S; R; n; angmean*180/pi; k; angmedian(ang)*180/pi]);

    fprintf('%8s %8.3f %8.3f %8.3f %8d %8.1f %8.2f %8.1f\n', 'Joint',Call,Sall,...
        Rall, nall, angmeanall*180/pi,kall, angmedian(ang(:))*180/pi);
    fprintf('\n');
end;

switch lower(method),
    case {'ww','watson-williams'},
        % From Mardia 2000, p.135
        if (~quiet && any(k < 1)),
            warning('anganova:LowConcentration',['Approximation is poor for kappa < 1.' ...
                ' Use likelihood ratio test']);
        end;

        if (~quiet),
            fprintf('High-concentration F test\n\nANOVA table\n');
            fprintf('%10s %5s %8s %8s %8s %8s %10s %10s\n', 'Source','d.f.','SS',...
                'MS','F', 'P','Modif F','Modif P');
        end;

        dfb = ngp-1;
        SSb = sum(R)-Rall;
        MSb = SSb/dfb;

        dfw = nall-ngp;
        SSw = nall - sum(R);
        MSw = SSw/dfw;

        F = MSb/MSw;
        Fmod = F*(1-3/8/kall);
        P = 1 - fcdf(F,dfb,dfw);
        Pmod = 1 - fcdf(Fmod,dfb,dfw);

        if (~quiet),
            fprintf('%10s %5d %8.3f %8.3f %8.3f %8.5f %8.3f %8.5f\n', 'Between',...
                dfb,SSb,MSb,F,P,Fmod,Pmod);
            fprintf('%10s %5d %8.3f %8.3f\n', 'Within',dfw,SSw,MSw);
        end;

    case {'lr','likelihood','likelihood ratio'},
        % From Mardia 2000, p.136
        angmean = atan2(S,C);
        angmeanall = atan2(Sall,Call);

        w = 2*kall*sum(R.*(1-cos(angmean-angmeanall)));

        U = 2/nall*(sum(R)^2 - Rall^2);

        c = 1/(1 - 1/8*kall^2 + ngp/(2*nall*kall^2));

        xsq = c*U;
        P = 1 - chi2cdf(xsq,ngp-1);

        if (~quiet),
            fprintf('Likelihood ratio test\n');
            fprintf('%8s %8s %8s %8s %8s %8s\n', 'w','U','c','Xsq','df','P');
            fprintf('%8.3f %8.3f %8.3f %8.3f %8d %8.5f\n', w,U,c,xsq,ngp-1,P);
        end;

    case {'concentration','conc'},
        % From Fisher 1993 p.131
        angmean = atan2(S,C);
        d = abs(sin(ang-repmat(angmean,[size(ang,1) 1])));

        dwithin = nansum(d)./n;
        dbetween = nansum(d(:))./nall;

        dfb = ngp-1;
        SSb = nansum(n.*(dwithin-dbetween).^2);
        MSb = SSb/dfb;

        dfw = nall-ngp;
        SSw = nansum(nansum((d - repmat(dwithin,[size(d,1) 1])).^2));
        MSw = SSw/dfw;

        F = MSb/MSw;
        P = 1 - fcdf(F,ngp-1,nall-ngp);

        if (~quiet),
            fprintf('Fisher test for equality of concentration\n\nANOVA table\n');
            fprintf('%10s %5s %8s %8s %8s %8s\n', 'Source','d.f.','SS',...
                'MS','F', 'P');
            fprintf('%10s %5d %8.3g %8.3g %8.3g %8.5g\n', 'Between',...
                dfb,SSb,MSb,F,P);
            fprintf('%10s %5d %8.3g %8.3g\n', 'Within',dfw,SSw,MSw);
        end;

    case {'npdist','mardia-watson-wheeler','mww'},
        % from Mardia 2000, p.157
        % also Batschelet 1981, p.101-105

        %we want to do something like the reverse of sorting - assign each
        %element a rank, based on where it would be if the angles were sorted
        %first sort them
        [rankval,rankinv] = sort(mod(ang(:),2*pi));
        %now ord is the rank of elements in the sorted list (just 1 to the
        %number of items, by definition, since they're sorted)
        ord = 1:numel(ang);
        %and we reverse the process by assigning the correct order value to
        %each element from the original matrix (with indices specified by
        %rankinv)
        rank(rankinv) = ord;

        rank = reshape(rank,size(ang));

        %find ties
        %based on how the sort function works, elements with the same value
        %will have the same rank
        tie = find(diff(rankval) == 0);
        if (~isempty(tie)),
            i = 1;
            nties = 0;
            while (i <= length(tie)),
                len = 1;
                while ((tie(i)+len <= n) && ...
                        (rankval(tie(i)+len) == rankval(tie(i)))),
                    len = len+1;
                end;

                % occasionally we get a weird rounding error where the
                % difference is zero,
                % but the values aren't directly equal.  Check for it here
                % so we don't get an infinite loop
                if (len == 1),
                    i = i+1;
                    continue;
                end;

                %shuffle the ranks in the ties randomly
                shuf0 = 1:len;
                shuffle = zeros(1,len);
                for j = 1:len-1,
                    mv = round(rand(1)*(len-j))+1;
                    shuffle(j) = shuf0(mv);
                    shuf0 = shuf0([1:mv-1 mv+1:end]);
                end;
                shuffle(len) = shuf0;

                %increase the ranks of items after this tie by the number of
                %elements in the tie
                rank(rankinv(tie(i):tie(i)+len-1)) = ...
                    tie(i) + shuffle-1;

                i = i+len-1;
                nties = nties + len;
            end;
            if (~quiet),
                warning('anganova:TieBroken','Broke %d ties.  Output will vary.',nties);
            end;
        end;

        %construct a new set of angles, where all the points are distributed
        %evenly around the circle, according to their ranks.  If points in
        %different groups are *still* separated, after we mangled them so
        %severely, we can conclude that they're from different populations
        beta = 2*pi*rank/nall;

        Rsq = sum(cos(beta)).^2 + sum(sin(beta)).^2;
        W = 2*sum(Rsq./n);                   % test statistic

        %approximately chi-squared distribution
        P = 1 - chi2cdf(W,2*(ngp-1));

        if (~quiet),
            fprintf('Mardia-Watson-Wheeler test:\n');
            fprintf('%8s %8s %8s\n','W','df','P');
            fprintf('%8.3f %8d %8.5f\n',W,2*(ngp-1),P);
        end;

    case {'median'},
        % from Fisher 1993, p.114
        med = angmedian(ang(:));

        m = sum(mod((ang - med)+pi,2*pi) < pi);
        M = sum(m);

        Pr = nall^2/(M*(nall-M))*sum(m.^2./n) - nall*M/(nall-M);

        P = 1 - chi2cdf(Pr,ngp-1);

        if (~quiet),
            fprintf('Fisher common median test\n');
            fprintf('%8s %8s %8s\n','Pr','df','P');
            fprintf('%8.3g %8d %8.5g\n',Pr,ngp-1,P);
        end;


    case {'npmean','watson'},
        %from Fisher 1993 p. 115-117
        if (~quiet && any(n < 25)),
            warning('anganova:SmallSampleSize','Test statistic is probably not good for n < 25.');
        end;

        %calculate the second moments
        rho2 = 1./n .* nansum(cos(2*(ang - repmat(angmean,[size(ang,1) 1]))));
        %and the dispersion (Fisher, 1993, eq 2.28)
        dsp = (1-rho2)./(2*(R./n).^2);

        %two tests, depending on whether the dispersions are comparable
        if (max(dsp)/min(dsp) <= 4),
            %dispersions are relatively similar
            Cp = nansum(n.*cos(angmean));
            Sp = nansum(n.*sin(angmean));

            R = sqrt(Cp^2 + Sp^2);

            dsp0 = nansum(n.*dsp)/nall;

            Y = 2*(nall - R)/dsp0;
        else
            %dispersions are different
            sigma2 = dsp./n;

            Cm = nansum(cos(angmean)./sigma2);
            Sm = nansum(sin(angmean)./sigma2);
            R = sqrt(Cm^2 + Sm^2);

            Y = 2*(nansum(1./sigma2) - R);
        end;

        P = 1 - chi2cdf(Y,ngp-1);

        if (~quiet),
            fprintf('Watson nonparametric means test\n');
            fprintf('%8s %8s %8s %8s\n','R','Y','df','P');
            fprintf('%8.3g %8.3g %8.1d %8.5g\n', R,Y,ngp-1,P);
        end;
end;

if (~quiet && (ngp > 2) && strcmp(display,'long')),
    %post-hoc tests
    %calculate number of tests
    ntests = sum(1:ngp-1);
    %Bonferroni reduction in alpha
    alphabonf = alpha/ntests;
    sigdif = cell(ngp,ngp+1);
    [sigdif{:}] = deal(' ');

    fprintf('\nPost-hoc test.  Bonferroni adjusted alpha = %f\n', alphabonf);
    fprintf('Unadjusted P values:\n');
    
    Pdif = NaN(ngp,ngp);
    for gp1 = 1:ngp,
        sigdif{gp1,1} = gp(gp1);
        fprintf(gptplt,gp(gp1));
        for gp2 = 1:gp1-1,
            P1 = anganova(ang(:,[gp1 gp2]),method, 'display','off');
            Pdif(gp1,gp2) = P1;
            if (P1 < alphabonf/10),
                sigdif{gp1,gp2+1} = '***';
            elseif (P1 < alphabonf/5),
                sigdif{gp1,gp2+1} = '**';
            elseif (P1 < alphabonf),
                sigdif{gp1,gp2+1} = '*';
            elseif (P1 < alpha),
                sigdif{gp1,gp2+1} = '-';
            end;
            fprintf(' %8.4g',Pdif(gp1,gp2)');
        end;
        fprintf('\n');
    end;
    fprintf(['%8s ' repmat([gptplt ' '],[1 ngp]) '\n'], ' ',gp);
    
    fprintf('Significant differences:\n');
    tplt = [gptplt repmat(' %8s',[1 ngp]) '\n'];
    sigdif = sigdif';
    fprintf(tplt, sigdif{:});
    tplt = repmat([gptplt ' '],[1 ngp]);
    tplt = ['%8s ' tplt(1:end-1) '\n'];
    fprintf(tplt,' ',gp);

    fprintf('*** = P < %f, ** = P < %f, * = P < %f, - = P < %f\n', ...
         alphabonf/10, alphabonf/5, alphabonf, alpha);
end;

if (nargout == 1),
    varargout = {P};
end;
