function [stats] = mediation_power(varargin)
% function [stats] = mediation([stats],[plots])
%
% Tor Wager, March 2006
%
% X, Y, M can be
% 1) vectors of observations
% 2) matrices of N columns (observations x N subjects)
% 3) cell arrays of length N (each cell is vector of obs. for one
% subject)
%
%
% columns of paths:
% 1 a   X -> M relationship
% 2 b   M -> Y relationship
% 3 cp  unmediated X -> Y relationship (residual)
% 4 c   X -> Y relationship
% 5 ab  mediated X -> Y by M (a * b)
%
% Examples:
%stats = mediation_power('reps',2,'N',20,'iter',20);
%
% stats = mediation_power('reps',23,'N',36,'iter',20,'p',.2,'alph',.001,'groupstd',.5); stats
% stats = mediation_power('multi','reps',reps,'N',60,'iter',10,'p',.4,'alph',.001,'groupstd',.33);

% --------------------------------------------
% Setup
% --------------------------------------------

% defaults
p = .4;  % proportion of variance transmitted x->m and m->y
N = 60;
iter = 50;
alph = .05;
reps = 1;       % replications
groupstd = .33;

domulti = 0; doplots = 0; dorobust = 1; verb = 1;
vnames = {'X' 'Y' 'M'};

for i = 1:length(varargin)
    if iscell(varargin{i})
        % do nothing; we have a non-string input
    else
        switch varargin{i}
            case 'multi', domulti = 1;
            case 'plots', doplots = 1;
            case 'robust', dorobust = 1;
            case 'verbose', verb = 1;
            case 'noverbose',verb = 0;
            case 'names', vnames = varargin{i+1};
                
            case 'p', p = varargin{i+1};
            case 'N', N = varargin{i+1};
            case 'iter', iter = varargin{i+1};  
            case 'reps', reps = varargin{i+1}; 
            case 'alph', alph = varargin{i+1};
            case 'groupstd', groupstd = varargin{i+1}; 
            otherwise
        end
    end
end


% multi mode: search over N
if domulti
    stats = [];
    reprep = 10:5:40
    for repno = 1:length(reprep)
        tmpstats = mediation_power('reps',reprep(repno),'N',N,'iter',iter,'p',p,'alph',alph,'groupstd',groupstd);
        stats.ind_runs(repno) = tmpstats;
        stats.all_power(repno,:) = tmpstats.power;
        disp(num2str(tmpstats.power))
    end
    
    figure('Color','w'); set(gca,'FontSize',16);
    plot(reprep,stats.all_power(:,5),'k-','LineWidth',2); 
    title('Power to detect mediation');
    xlabel('Subjects'); ylabel('Power');
    return
end



obs = N .* reps;    % total observations to create, and reshape if necessary
if reps > 1, statstring = 'stats';, else, statstring = 'boot1';, end

if verb
    fprintf(1,'Mediation power\n\np: %3.2f, Iterations: %3.0f\n',p,iter);
    fprintf(1,'---------------------------------\n')
end
paths = zeros(iter,5); pp = paths;
verbstring = 'verbose';

% --------------------------------------------
% Generate data
% --------------------------------------------
% for i = 1:100
%     xym = mvnrnd([0 0 0],[1 .5 .5; .5 1 .5; .5 .5 1],1000);
%     [paths(i,:)] = mediation(xym(:,1),xym(:,2),xym(:,3));
% end
% 
% % x,y unrelated; both x and y contribute to m
% % ends up looking like suppression
% for i = 1:100
%     xy = mvnrnd([0 0],[1 0; 0 1],1000);
%     m = xy(:,1)+normrnd(0,1,1000,1) + xy(:,2)+normrnd(0,1,1000,1);
%     [paths(i,:)] = mediation(xy(:,1),xy(:,2),zscore(m));
% end
% 
% for i = 1:100
%     x = mvnrnd(0,1,1000);
%     xym = mvnrnd([0 0 0],[1 .5 .5; .5 1 .5; .5 .5 1],1000);
%     [paths(i,:)] = mediation(xym(:,1),xym(:,2),xym(:,3));
% end


% all mediated, no direct effect.  single level

if reps > 1, gx = (1:N)'; gx = zscore(gx); , end
p2 = p .^ 2;

for i = 1:iter

    xym = mvnrnd([0 0 0],[1 p2 p; p2 1 p; p p 1],obs);
    x = xym(:,1); y = xym(:,2); m = xym(:,3);
    

    docheck = 0;

    % reshape for group analysis if necessary
    if reps > 1
        x = reshape(x,N,reps);
        y = reshape(y,N,reps);
        m = reshape(m,N,reps);

        if docheck
            clear bb
            for ii=1:reps
                bb(:,ii) = pinv([x(:,ii) ones(size(x,1),1)]) * y(:,ii);
            end
            
            bb = bb(1,:);bb1 = std(bb); fprintf(1,'\nBefore: std(b) = %3.2f     ',std(bb));
        end

        % add between-subjects noise
        n = repmat(normrnd(0,groupstd,1,reps),N,1);
        m = m + x .* n;
        n = repmat(normrnd(0,groupstd,1,reps),N,1);
        y = y + m .* n;

        if docheck
            clear bb
            for ii=1:reps
                bb(:,ii) = pinv([x(:,ii) ones(size(x,1),1)]) * y(:,ii);
            end
            bb = bb(1,:); fprintf(1,'After: std(b) = %3.2f    ',std(bb));
            fprintf(1,'Diff: %3.2f\n    ',std(bb) - bb1);
        end

    end
    
    %cors(:,:,i) = corrcoef([x m y]);
    
    [mypath,stats2] = mediation(x,y,m,statstring,verbstring);
    
    if i == 1
        fprintf(1,'---------------------------------\n')
        if verb, fprintf(1,'Simulating power in %3.0f iterations: %03d',iter,1);, end
    else
        if verb, fprintf(1,'\b\b\b%03d',i);,end
    end
    
    if reps > 1, paths(i,:) = stats2.mean;, else, paths(i,:) = mypath;, end
    pp(i,:) = stats2.p;
    
    verbstring = 'xxx'; % turn off verbose
    
end

% --------------------------------------------
% Summary
% --------------------------------------------

sig = pp < alph;  %sign(paths) .* (pp < alph);  % p-values are already 2-tailed
pow = sum(sig) ./ length(sig);
%cor = squeeze(mean(cors,3));

stats = struct('power',pow,'prop',p,'obs',N,'iter',iter, 'replications',reps,...
    'alpha',alph,'paths',paths,'sig',sig,'p',pp);

stats.names = stats2.names;



return











% --------------------------------------------
%
%
% Sub-functions
%
%
% --------------------------------------------

function stats = getstats(bp)

%stats.names = {'a' 'b' 'c''' 'c' 'ab'};
%stats.mean = mean(bp);
%stats.ste = std(bp);
%stats.ci95 = [prctile(bp,97.5); prctile(bp,2.5)];
%stats.ci90 = [prctile(bp,95); prctile(bp,5)];
%stats.ci99 = [prctile(bp,99.5); prctile(bp,.5)];
stats.p = 2.*(min(sum(bp<=0),sum(bp>=0)) ./ size(bp,1));

return

function plot_hists(bootpaths,vnames);

[h,xx] = hist(bootpaths(:),50);

myperc = mean(xx).*.1;
xlimit = [min(xx)-myperc max(xx)+myperc];
a=bootpaths(:,1); b=bootpaths(:,2);cp=bootpaths(:,3);c=bootpaths(:,4);ab=bootpaths(:,5);
tor_fig(1,5);
shaded_hist(a,xx);
title(['a: ' vnames{1} '->' vnames{3}]); set(gca,'XLim',xlimit);

subplot(1,5,2);shaded_hist(b,xx);
title(['b: ' vnames{3} '->' vnames{2}]);set(gca,'XLim',xlimit);

subplot(1,5,3);shaded_hist(c,xx);
title(['c: ' vnames{1} '->' vnames{2}]);set(gca,'XLim',xlimit);

subplot(1,5,4);shaded_hist(cp,xx);
title(['c'':' vnames{1} '->' vnames{2}]);set(gca,'XLim',xlimit);

subplot(1,5,5);shaded_hist(ab,xx);
title(['ab: ' vnames{1} '->' vnames{2}]);set(gca,'XLim',xlimit);
drawnow

return




function shaded_hist(a,xx)

h=hist(a,xx); han = bar(xx,h); set(han,'FaceColor',[.7 .7 .7]); set(han,'EdgeColor',[.7 .7 .7]);
wh = (xx > prctile(a,97.5) | xx < prctile(a,2.5));
h(~wh) = 0;
if any(wh'), hold on; han = bar(xx,h); set(han,'FaceColor',[.3 .3 .3],'EdgeColor',[.3 .3 .3]);, end

return



function plot_slopes(abetas,bbetas,cbetas,cpbetas,X,M,vnames)
% abetas and others are matrices: rows = [slope, intercept], cols =
% subjects


tor_fig(1,4);
plot_slopes_subfunction(abetas,X);  % draws a line for each subject.
xlabel(vnames{1}), ylabel(vnames{3}); title('a: X->M');

subplot(1,4,2);
plot_slopes_subfunction(bbetas,M);  % draws a line for each subject.
xlabel(vnames{3}), ylabel(vnames{2}); title('b: M->Y controlling X');

subplot(1,4,3);
plot_slopes_subfunction(cbetas,X);  % draws a line for each subject.
xlabel(vnames{1}), ylabel(vnames{2}); title('c: X->Y');

subplot(1,4,4);
plot_slopes_subfunction(cpbetas,X);  % draws a line for each subject.
xlabel(vnames{1}), ylabel(vnames{2}); title('c'': X->Y controlling M');

return


function plot_slopes_subfunction(betas,X)
    for i = 1:size(betas,2)    % for each subject
        if iscell(X), x = X{i};, else, x = X(:,i);, end

        % 95% conf. interval for x
        mx = mean(x); s=std(x)*tinv(.975,length(x)-1); x = [mx - s mx + s];
        y = betas(2,i) + betas(1,i) * x;
        plot(x,y,'k');
    end
    
    drawnow
return


    

