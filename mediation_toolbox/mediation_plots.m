% Plotting function for mediation output.
% Several kinds of plots can be created.
%
% This is used in mediation.m
%
% See also mediation_scatterplots.
% Tor Wager, Jan 2009
%
% Usage:
% -------------------------------------------------------------------------
% Plot individual slopes of regressions, using conf. interval for x range
% -------------------------------------------------------------------------
%
% mediation_plots(stats2, 'slopes', 'noint')
% mediation_plots(stats2, 'slopes', 'noint', 'color', [.3 .3 1])
%
%
% -------------------------------------------------------------------------
% Plot the covariance between a and b paths in a multi-level model
% -------------------------------------------------------------------------
% mediation_plots(stats2, 'abcov')

function mediation_plots(stats2, method_name, varargin)

doslopes = 1;
categ = [0 0 0];  % 1 for categorical, 0 otherwise
autocat = 1;          % automatically choose categorical vs. continuous

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % General Defaults
            case {'noslopes'}, doslopes = 0;
                 
            % Other arguments used in subfunctions
            case {'noint', 'nointercept'}
            case 'noind', 
            case 'samefig',             
            case {'weight', 'weighted', 'var', 's2'}
            case 'noweight'
            case 'color'
                
            otherwise
                disp('Warning! Unknown input string option.');
        end
    end
end

if autocat
    indx = 1;
    for vars = {'X' 'Y' 'M'}
        myvar = stats2.inputOptions.(vars{1});
        len = length(unique(cat(1, myvar{:})));

        if len < 3, categ(indx) = 1; end
    end
end
    
switch method_name
    case {'slope', 'slopes', 'line', 'lines'}
        slope_plots(stats2, varargin{:})

    case 'abcov'
        ab_cov_plot(stats2, varargin{:});
        
    otherwise error('Unknown method name for mediation_plots.m');
end

end




function slope_plots(stats2, varargin)

newfig = 1;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % General Defaults
            case 'samefig', newfig = 0;

            otherwise
                disp('Warning! Unknown input string option.');

        end
    end
end

% abetas and others are matrices: rows = [slope, intercept], cols =
% subjects

N = stats2.inputOptions.N;
npaths = stats2.inputOptions.npaths;

% abetas, etc.: each column represents a regression eq. (for a paths) or a mediator (for
% b paths)
abetas = stats2.abetas;
bbetas = stats2.bbetas;
cbetas = stats2.cbetas;
cpbetas = stats2.cpbetas;

X = stats2.inputOptions.X;
M = stats2.inputOptions.M;
additionalM = stats2.inputOptions.additionalM;

if iscell(additionalM), num_additionalM = size(additionalM{1}, 2); 
else num_additionalM = size(additionalM, 2); 
end

vnames = stats2.inputOptions.vnames;

if iscell(abetas)
    whgood = true(N, 1);
    for i = 1:length(abetas)
        if isempty(abetas{i}) || isempty(bbetas{i}) || isempty(cbetas{i}) || isempty(cpbetas{i})
            whgood(i) = 0;
        end
    end

    abetas = cat(2, abetas{:});
    bbetas = cat(2, bbetas{:});
    
    % multiple mediators, reformat
    for i = 1:num_additionalM + 1
       tmp_abetas{i} = abetas(:, i : num_additionalM + 1 : end); 
       tmp_bbetas{i} = bbetas([i end], :); 
    end
    abetas = tmp_abetas;
    bbetas = tmp_bbetas;
    
    cbetas = cat(2, cbetas{:});
    cpbetas = cat(2, cpbetas{:});
else
    % need to extend for multiple mediators ***
    whgood = ~any(isnan([abetas' bbetas' cbetas' cpbetas']), 2);
    
    abetas = {abetas}; % for consistency with multi-level/multi-mediator
    bbetas = {bbetas};
end

% stats2 is varargin{1}
% w is [a b c' c ab] weights, N rows
w = ones(N, npaths);

% we have stats structure with (maybe) weights
if isfield(stats2, 'w')
    w = stats2.w;
    w = w .* N;
end

w = w(whgood,:);

ncols = 4;
nrows = num_additionalM + 1;

fh = create_figure('Slope_Plot', nrows, ncols);

for wh_row = 1:nrows % each row is a different first-level mediator
    
    prevrows = ncols * (wh_row - 1);
    
    % get the correct weights and names for each mediator
    if wh_row == 1
        myw = w(:, 1:2);
        wh_mediator_name = 3;
        
    else
        % This is this way because of the specific ordering of columns
        % containing stats for additional mediators.
        myw = w(:, 5 + 3 * (wh_row-1) - 2 : 5 + 3 * (wh_row-1) - 1); 
        wh_mediator_name = wh_row - 1 + 3;
    end
    
    
subplot(nrows, ncols, 1 + prevrows); hold on;
plot_individual_slopes(abetas{wh_row}, X, myw(:,1), varargin{:});
xlabel(vnames{1}), ylabel(vnames{wh_mediator_name});
if wh_row == 1, title('a: X->M'); end

subplot(nrows, ncols, 2 + prevrows); hold on;
plot_individual_slopes(bbetas{wh_row}, M, myw(:,2), varargin{:});
xlabel(vnames{wh_mediator_name}), ylabel(vnames{2});
if wh_row == 1, title('b: M->Y controlling X'); end

if wh_row == 1  % these are only defined for the first row
    subplot(nrows, ncols, 3); hold on;
    plot_individual_slopes(cpbetas, X, w(:,3), varargin{:});
    xlabel(vnames{1}), ylabel(vnames{2});
    title('c'': X->Y controlling M');

    subplot(nrows, ncols, 4); hold on;
    plot_individual_slopes(cbetas, X, w(:,4), varargin{:});
    xlabel(vnames{1}), ylabel(vnames{2});
    title('c: X->Y');
else
    subplot(nrows, ncols, 3 + prevrows); axis off
    subplot(nrows, ncols, 4 + prevrows); axis off
end

end % rows (mediators)


end


function plot_individual_slopes(betas, X, w, varargin)

doslopes = 1;

nointercept = 0; % suppress intercept
newfig = 1;
weight_option = 'weighted';
doind = 1;

baseindlinecolor = [0 0 1];

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % General Defaults
            case {'noint', 'nointercept'}, nointercept = 1;
            case 'noind', doind = 0;
            case 'samefig', newfig = 0;

                % Estimation Defaults
            case {'weight', 'weighted', 'var', 's2'}, weight_option = 'weighted';
            case 'noweight', weight_option = 'unweighted';

            case 'color', baseindlinecolor = varargin{i + 1};

            otherwise
                disp('Warning! Unknown input string option.');

        end
    end
end

% sort so that we plot from lowest to highest weights
[w, sorti] = sort(w);
betas = betas(:,sorti);
if iscell(X), X = X(sorti); else X = X(:,sorti); end

% get colors based on weights
N = size(betas, 2);
minwt = .4;     % make sure all lines are visible.
w = rescale_range(w, [minwt 1]);

% line widths: median split, top half gets 2.
linew = (w(:,1) >= median(w(:,1))) + 1;

%w = 1 - w;  % for colors, 0 is black
white = [1 1 1];
colors = (repmat(w, 1, 3) .* repmat(baseindlinecolor, N, 1)) + repmat((1 - w), 1, 3) .* repmat( white, N, 1 );
colors(colors < 0) = 0; % include to remove rounding error

% do not plot intercepts, if asked for
if nointercept, betas(2, :) = nanmean(betas(2, :), 2); end

%all_mx = zeros(N, 1);
allx = [];
for i = 1:N
    if iscell(X), x = X{i}; else x = X(:,i); end

    % 95% conf. interval for x
    [nanvec, x] = nanremove(x);
    mx = mean(x);

    %all_mx(i) = mx;
    s = std(x) * tinv(.975, length(x)-1);
    x = [mx - s mx + s];
    y = betas(2, i) + betas(1, i) * x;
    plot(x, y, '-', 'Color', colors(i,:), 'LineWidth', linew(i));

    allx = [allx x];
end

% set x range
xminval = min(allx) - .05 * range(allx);
xmaxval = max(allx) + .05 * range(allx);
set(gca, 'XLim', [xminval xmaxval]);
plot_mean_line(betas, w);


end


function rx = rescale_range(x, y)
% re-scale x to range of y
m = range(y)./range(x);

if isinf(m)
    % no range/do not rescale
    rx = x;
else
    b = y(1) - m * x(1);
    rx = m*x + b;
end
end


function plot_mean_line(betas, w)
% mean line
wmean = @(paths, w) ((w ./ sum(w))'*paths)';
means = wmean(betas', w); % should be slope then intercept
x = get(gca, 'XLim');
y = means(2) + means(1) * x;

% SE lines (bootstrapped; multilevel)
xx = linspace(x(1), x(2), 50);

boots = 1000;
means =  bootstrp(boots, wmean, betas', w);
yy = zeros(boots, 50);
for i = 1:boots
    yy(i,:) = means(i, 1) * xx + means(i, 2);
end

ub = prctile(yy, 95);
plot(xx, ub, 'k', 'LineWidth', 2);
lb = prctile(yy, 5);
plot(xx, lb, 'k', 'LineWidth', 2);

fhan = fill([xx xx(end:-1:1)], [lb ub(end:-1:1)], [.5 .5 .5]);
set(fhan, 'FaceAlpha', .3)

plot(x, y, '-', 'Color', 'k', 'LineWidth', 3);

% xb = mean(x); s = resid. std. stat.se(2) is se of slope
%se_mean = (stat.se(2).^2 * (xx-xb).^2 + s.^2/length(x)) .^ .5;

drawnow

end





function ab_cov_plot(stats2, varargin)

    if ~isfield(stats2, 'analysisname'), disp('Cannot find stats2.analysisname. Skipping a,b cov plot.'); return, end
    
        
    if ~strcmp(stats2.analysisname, 'Multi-level model')
        disp('Model does not appeart to be a multi-level model.  Skipping a,b cov plot.');
        return
    end
    
    if ~isfield(stats2, 'ebayes_bstar')
        disp('Cannot find ebayes_bstar field in stats output, so cannot create cov(a,b) scatterplot');
        return
    end
    
    create_figure('cov(a,b)'); 
    
    
    plot_correlation_samefig(stats2.ebayes_bstar(:, 1), stats2.ebayes_bstar(:, 2));
    set(gca, 'FontSize', 24);
    
    xlabel('Path a (Empirical Bayes)');
    ylabel('Path b (Empirical Bayes)');

end

