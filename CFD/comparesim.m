function comparesim(files,varargin)
% function comparesim(files,varargin)
% Utility function for comparing various simulations.  Takes a cell array
% of simulation file names.

opt.nphase = 2;             % number of different phase values to show
opt.col = 'bgrcmk';         % colors for the different simulations
opt.stridelength = false;   % Show time normalized by the activation period, or not.
opt.spaceamp = 0;           % Space amplitude plots out vertically
opt.bodylen = 4*pi;         % Body length
opt.sidebyside = false;     % Show plots on top of each other, or side-by-side
opt.labels = {};            % labels for each of the data files

opt = parsevarargin(opt,varargin,2,'typecheck',false);

N = length(files);

if (ischar(opt.col)),
    opt.col = opt.col(:);
end;

s = cell(1,N);
t = cell(1,N);
comspeed = cell(1,N);
comx = cell(1,N);
comy = cell(1,N);
xm = cell(1,N);
ym = cell(1,N);
ampcont = cell(1,N);
env = cell(1,N);
Ithrust = cell(1,N);
Idrag = cell(1,N);
worktot = cell(1,N);
skinfricthrust = cell(1,N);
skinfricdrag = cell(1,N);
rmslateralt = cell(1,N);
rmslateraln = cell(1,N);
freq = zeros(1,N);
for i = 1:N,
    F = load(files{i}, 't','freq','xm','ym','xl','yl','xr','yr','xn','yn',...
        'um','vm','ul','vl','ur','vr','un','vn','good');
    [pn,fn] = fileparts(files{i});
    
    F = joinstructfields(F,load(fullfile(pn,[fn '_analysis.mat']), 'ampcont','comspeed', 'comx','comy',...
        'worktot','workpos','workneg',...
        'workposact','worknegact','Stress','s','steadycycle','cyclenum', 'ps', ...
        'swimvecx','swimvecy'));
    
    if (isfield(F,'good')),
        F = resizeallvars(2,F.good, 'struct',F);
    end;
    
    freq(i) = F.freq;
    t{i} = F.t;
    s{i} = F.s;
    comspeed{i} = F.comspeed;
    comx{i} = F.comx;
    comy{i} = F.comy;
    xm{i} = F.xm;
    ym{i} = F.ym;
    
    env{i} = F.ampcont(:,F.cyclenum >= F.steadycycle);
    
    c = max(F.cyclenum) - 1;
    a = first(F.cyclenum == c);
    b = last(F.cyclenum == c) + 1;
    ind = round(linspace(a,b, opt.nphase+1));
    ampcont{i} = F.ampcont(:,ind(1:end-1));
    
    worktot{i} = F.worktot(:,F.steadycycle:end)*F.ps;
    
    if (isfield(F,'Stress'))
        F.Stress = fluidstressnearboundary(makeib(F),{},F.swimvecx,F.swimvecy, 'fluidvals',F.Stress);
        F.Stress = stress2thrustdrag(F.t,F.cyclenum,F.Stress);
    
        Ithrust{i} = sum(F.Stress.Ithrustt(:,F.steadycycle:end,:) + ...
            F.Stress.Ithrustn(:,F.steadycycle:end,:),3);
        Idrag{i} = sum(F.Stress.Idragt(:,F.steadycycle:end,:) + ...
            F.Stress.Idragn(:,F.steadycycle:end,:),3);
        skinfricthrust{i} = sum(F.Stress.Ithrustt(:,F.steadycycle:end,:),3);
        skinfricdrag{i} = sum(F.Stress.Idragt(:,F.steadycycle:end,:),3);
        
        rmslateralt1 = NaN(size(F.Stress.Ithrustt,1),size(F.Stress.Ithrustt,2));
        rmslateraln1 = NaN(size(F.Stress.Ithrustt,1),size(F.Stress.Ithrustt,2));
        for c = F.steadycycle:max(F.cyclenum)-1,
            iscycle = F.cyclenum == c;
            rmslateralt1(:,c) = sqrt(nanmean(sum(F.Stress.Flateralt(:,iscycle,:).^2,3),2));
            rmslateraln1(:,c) = sqrt(nanmean(sum(F.Stress.Flateraln(:,iscycle,:).^2,3),2));
        end;
        
        rmslateralt{i} = rmslateralt1;
        rmslateraln{i} = rmslateraln1;
    end;
end;

s = catuneven(3,s{:});
t = catuneven(3,t{:});
comspeed = catuneven(3,comspeed{:});
comx = catuneven(3,comx{:});
comy = catuneven(3,comy{:});
xm = catuneven(3,xm{:});
ym = catuneven(3,ym{:});
ampcont = catuneven(3,ampcont{:});
env = catuneven(3,env{:});
Ithrust = catuneven(3,Ithrust{:});
Idrag = catuneven(3,Idrag{:});
worktot = catuneven(3,worktot{:});
skinfricthrust = catuneven(3,skinfricthrust{:});
skinfricdrag = catuneven(3,skinfricdrag{:});
rmslateralt = catuneven(3,rmslateralt{:});
rmslateraln = catuneven(3,rmslateraln{:});

figure(1);
clf;

if (opt.sidebyside),
    h(1,1) = subplot(2,1,1);
    h(1,2:N) = h(1,1);
    for j = 1:N,
        h(2,j) = subplot(2,N,N+j);
    end;
else
    h(1,1) = subplot(3,1,1);
    h(2,1) = subplot(3,1,2:3);
    h(1:2,2:N) = repmat(h(1:2,1),[1 N-1]);
end;

figure(2);
clf;
h(3,1) = subplot(3,1,1);
h(4,1) = subplot(3,1,2);
h(5,1) = subplot(3,1,3);
h(3:5,2:N) = repmat(h(3:5,1),[1 N-1]);

pt = linspace(1,size(s,1), size(Ithrust,1)+1);
pt = (pt(1:end-1) + pt(2:end))/2;

for i = 1:N,
    if (opt.stridelength),
        addplot(h(1,i), t(:,:,i)*freq(i),comspeed(:,:,i)/freq(i)/opt.bodylen, ...
        '-','LineWidth',2,'Color',opt.col(i,:));
    else
        addplot(h(1,i), t(:,:,i),comspeed(:,:,i)/opt.bodylen,'-','LineWidth',2,'Color',opt.col(i,:));
    end;
    
    addplot(h(2,i), s(:,:,i)/opt.bodylen,ampcont(:,:,i)/opt.bodylen+opt.spaceamp*i, '-', 'Color',opt.col(i,:));
    addplot(h(2,i), s(:,:,i)/opt.bodylen,max(env(:,:,i)/opt.bodylen+opt.spaceamp*i,[],2), '-', ...
        'LineWidth',2,'Color',opt.col(i,:));
    
    addplot(h(3,i), s(:,:,i)/opt.bodylen,nanmean(worktot(:,:,i),2), '-', ...
        'LineWidth',2,'Color',opt.col(i,:));
    
    addplot(h(4,i), s(pt,:,i)/opt.bodylen,nanmean(Ithrust(:,:,i),2), '-', ...
        'LineWidth',2,'Color',opt.col(i,:));
    addplot(h(4,i), s(pt,:,i)/opt.bodylen,nanmean(Idrag(:,:,i),2), '--', ...
        'LineWidth',2,'Color',opt.col(i,:));

    addplot(h(5,i), s(pt,:,i)/opt.bodylen,nanmean(rmslateraln(:,:,i)+rmslateralt(:,:,i),2), ...
        'LineWidth',2,'Color',opt.col(i,:));
end;

for i = 1:numel(h),
    axis(h(i),'tight');
end;
if (opt.sidebyside),
    ax1 = axis(h(1,1));
    ax2 = axis(h(1,2));
    ax1([1 3]) = min([ax1([1 3]); ax2([1 3])]);
    ax1([2 4]) = max([ax1([2 4]); ax2([2 4])]);
    
    axis(h(1,1),ax1);
    axis(h(1,2),ax1);

    ax1 = axis(h(2,1));
    ax2 = axis(h(2,2));
    ax1([1 3]) = min([ax1([1 3]); ax2([1 3])]);
    ax1([2 4]) = max([ax1([2 4]); ax2([2 4])]);
    
    axis(h(2,1),ax1);
    axis(h(2,2),ax1);
end;

if (opt.stridelength),
    ylabel(h(1),'COM speed (L/cycle)');
    xlabel(h(1),'Time (cycles)');
else
    ylabel(h(1,1),'Speed (L/s)');
    xlabel(h(1,1),'Time (s)');
end;
xlabel(h(2,1),'Position along body (L)');
ylabel(h(2,1),'Amplitude (L)');
xlabel(h(4),'Position along body (cm)');
ylabel(h(3),'Total muscle work');
ylabel(h(4),'Axial fluid impulse');
ylabel(h(5),'RMS lateral forces');

if (~isempty(opt.labels))
    legend(h(1,1),opt.labels,'Orientation','horizontal','Location','NorthOutside');
    legend(h(3),opt.labels,'Orientation','horizontal','Location','NorthOutside');
end;

fprintf('Skin friction thrust fraction: ');
fprintf('%g ',nanmean(flatten(skinfricthrust./Ithrust,1:2)));
fprintf('\n');
fprintf('Skin friction drag fraction: ');
fprintf('%g ',nanmean(flatten(skinfricdrag./Idrag,1:2)));
fprintf('\n');

fprintf('Skin friction lateral force fraction: ');
fprintf('%g ',nanmean(flatten(rmslateralt./(rmslateraln+rmslateralt),1:2)));
fprintf('\n');



    